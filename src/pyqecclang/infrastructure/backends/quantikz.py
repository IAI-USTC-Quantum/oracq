"""RIR 模块的 quantikz 线路导出：保持模块结构，不展开调用与重复。"""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence
from typing import cast

from pyqecclang.infrastructure.builder import Operation
from pyqecclang.infrastructure.ir import (
    Adjoint,
    Call,
    Control,
    Instruction,
    Load,
    Module,
    Primitive,
    Program,
    Ref,
    Repeat,
    Store,
    ValidationError,
)

_GATE_LABELS = {"h": "H", "x": "X", "y": "Y", "z": "Z", "s": "S", "t": "T"}
_ROTATION_LABELS = {"rx": "R_x", "ry": "R_y", "rz": "R_z", "phase": "P"}
_name_counter = itertools.count(1)


def _escape(name: str) -> str:
    """转义名字中的下划线，供 LaTeX 数学模式使用。"""
    return name.replace("_", r"\_")


def _angle_text(angle: float) -> str:
    """小整数倍 pi 排成分数线形式，其余按十进制。"""
    ratio = angle / math.pi
    for q in (1, 2, 3, 4, 6, 8, 12, 16):
        p = round(ratio * q)
        if abs(ratio * q - p) < 1e-9:
            if p == 0:
                return "0"
            sign = "-" if p < 0 else ""
            p = abs(p)
            if q == 1:
                return sign + ("\\pi" if p == 1 else f"{p}\\pi")
            num = "" if p == 1 else str(p)
            return sign + f"\\frac{{{num}\\pi}}{{{q}}}"
    return f"{angle:.6g}"


def _runs(rows: list[int]) -> list[list[int]]:
    """把行号列表切成若干连续段（保持给定顺序）。"""
    result: list[list[int]] = []
    for row in rows:
        if result and row == result[-1][-1] + 1:
            result[-1].append(row)
        else:
            result.append([row])
    return result


class _Renderer:
    """把入口模块逐节点摆放到 quantikz 矩阵列，并登记分组框。"""

    def __init__(self, program: Program, wire_labels: bool, name: str) -> None:
        """绑定程序、线标签选项与模块名，并初始化量子线表。"""
        self.program: Program = program
        self.module: Module = program.main
        self.wire_labels: bool = wire_labels
        self.name: str = name
        self.wires: list[tuple[str, int]] = []
        for register in (*self.module.registers, *self.module.locals):
            self.wires.extend((register.name, bit) for bit in range(register.type.width))
        self.row_of: dict[tuple[str, int], int] = {wire: row for row, wire in enumerate(self.wires)}
        self.columns: list[list[str]] = []
        self.groups: list[tuple[list[tuple[int, int]], str]] = []

    def rows(self, ref: Ref) -> list[int]:
        """把视图展开为其覆盖的矩阵行号列表。

        Args:
            ref: 待展开的寄存器视图。

        Returns:
            list[int]: 视图覆盖的量子线行号，按各段低位到高位的顺序排列。
        """
        result: list[int] = []
        for span in ref.parts:
            result.extend(
                self.row_of[(span.register, span.start + offset)] for offset in range(span.width)
            )
        return result

    def new_column(self) -> list[str]:
        """追加一列空单元格并返回该列。

        Returns:
            list[str]: 新追加的列；每个元素对应一条量子线的单元格。
        """
        self.columns.append([""] * len(self.wires))
        return self.columns[-1]

    def place_controls(
        self, column: list[str], target_rows: Sequence[int], controls: tuple[tuple[int, int], ...]
    ) -> None:
        """在列中放置控制点，指向各目标段的首行。

        Args:
            column: 待写入控制点标记的矩阵列。
            target_rows: 控制指向的目标行号集合，取其最小行作为指向端。
            controls: (行号, 生效位值) 元组序列；位值为 0 时画空心反控制点。
        """
        top = min(target_rows)
        for row, value in controls:
            mark = "\\ctrl" if value else "\\octrl"
            column[row] = f"{mark}{{{top - row}}}"

    def box(
        self,
        rows: list[int],
        label: str,
        controls: tuple[tuple[int, int], ...],
        *,
        dashed: bool = False,
    ) -> None:
        """把行号的每个连续段画成盒子门，按需添加控制点。

        Args:
            rows: 盒子覆盖的行号；连续段各自成盒并各占一列。
            label: 盒内的 LaTeX 标签。
            controls: 追加到每个盒子的控制点。
            dashed: 为真时盒子用虚线样式，表示开放声明等未绑定结构。
        """
        for run in _runs(rows):
            column = self.new_column()
            style = ",style={dashed}" if dashed else ""
            column[run[0]] = f"\\gate[wires={len(run)}{style}]{{{label}}}"
            if controls:
                self.place_controls(column, run, controls)

    def emit_primitive(self, node: Primitive, controls: tuple[tuple[int, int], ...]) -> None:
        """发射单个基元门：单比特门、CNOT、SWAP、加常数或全局相位。

        Args:
            node: 待绘制的基元指令节点。
            controls: 从外层 ``Control`` 节点继承的控制点。
        """
        op = node.op
        if op in _GATE_LABELS or op in _ROTATION_LABELS:
            label = _GATE_LABELS.get(op)
            if op in _ROTATION_LABELS:
                # 校验保证旋转类门的角度存在。
                label = f"{_ROTATION_LABELS[op]}({_angle_text(cast(float, node.angle))})"
            column = self.new_column()
            for row in self.rows(node.operands[0]):
                column[row] = f"\\gate{{{label}}}"
            if controls:
                self.place_controls(column, self.rows(node.operands[0]), controls)
            return
        if op == "xor":
            for source_row, target_row in zip(
                self.rows(node.operands[0]), self.rows(node.operands[1]), strict=True
            ):
                column = self.new_column()
                column[source_row] = f"\\ctrl{{{target_row - source_row}}}"
                column[target_row] = "\\targ{}"
                if controls:
                    self.place_controls(column, (target_row,), controls)
            return
        if op == "swap":
            for first, second in zip(
                self.rows(node.operands[0]), self.rows(node.operands[1]), strict=True
            ):
                column = self.new_column()
                column[first] = "\\swap"
                column[second] = "\\targX{}"
                if controls:
                    self.place_controls(column, (first, second), controls)
            return
        if op == "add_const":
            self.box(self.rows(node.operands[0]), f"+{node.value}", controls)
            return
        if op == "gphase":
            self.box(
                list(range(len(self.wires))),
                # 校验保证 gphase 的角度存在。
                f"e^{{i{_angle_text(cast(float, node.angle))}}}",
                controls,
            )
            return
        raise ValidationError(f"quantikz 暂不支持的门操作：{op}")

    def emit_load(self, node: Load, controls: tuple[tuple[int, int], ...]) -> None:
        """把 QRAM 查询画成资源盒子。

        Args:
            node: 待绘制的 QRAM 读取指令。
            controls: 从外层 ``Control`` 节点继承的控制点。
        """
        rows = self.rows(node.address) + self.rows(node.data)
        self.box(rows, f"\\mathrm{{QRAM}}_{{\\mathtt{{{_escape(node.resource)}}}}}", controls)

    def emit_store(self, node: Store, controls: tuple[tuple[int, int], ...]) -> None:
        """把 QRAM 随机写画成带上标 ``w`` 的资源盒子。

        Args:
            node: 待绘制的 QRAM 写入指令。
            controls: 从外层 ``Control`` 节点继承的控制点。
        """
        rows = self.rows(node.address) + self.rows(node.data)
        self.box(
            rows,
            f"\\mathrm{{QRAM}}^{{\\mathrm{{w}}}}_{{\\mathtt{{{_escape(node.resource)}}}}}",
            controls,
        )

    def emit_call(self, node: Call, controls: tuple[tuple[int, int], ...]) -> None:
        """把模块调用画成模块名盒子；开放声明用虚线。

        Args:
            node: 待绘制的模块调用指令。
            controls: 从外层 ``Control`` 节点继承的控制点。
        """
        rows: list[int] = []
        for argument in node.arguments:
            rows.extend(self.rows(argument))
        target = self.program.module_map.get(node.module)
        dashed = target is not None and target.body is None
        self.box(rows, f"\\mathtt{{{_escape(node.module)}}}", controls, dashed=dashed)

    def emit_grouped(
        self, node: Repeat | Adjoint, controls: tuple[tuple[int, int], ...], label: str
    ) -> None:
        """发射组内指令后，为其占用的单元登记一个分组标注框。

        Args:
            node: 组内待发射的 ``Repeat`` 或 ``Adjoint`` 结构节点。
            controls: 从外层 ``Control`` 节点继承的控制点。
            label: 分组框上方的 LaTeX 标注，如重复次数或 dagger。
        """
        start = len(self.columns)
        self.emit(node.body, controls)
        steps = len(self.columns) - start
        used = [
            (start + offset, row)
            for offset, column in enumerate(self.columns[start:])
            for row, cell in enumerate(column)
            if cell
        ]
        if not used:
            return
        if steps == 1:
            # 单盒子折叠为上标标记（DB^\dagger 形式更紧凑）。
            column = self.columns[start]
            filled = [(row, cell) for row, cell in enumerate(column) if cell]
            if len(filled) == 1 and filled[0][1].startswith("\\gate") and filled[0][1].endswith("}"):
                row, cell = filled[0]
                column[row] = cell[:-1] + f"^{{{label}}}}}"
                return
        self.groups.append((used, label))

    def emit(
        self, body: tuple[Instruction, ...], controls: tuple[tuple[int, int], ...] = ()
    ) -> None:
        """按节点类型分派发射；结构节点递归或折叠为分组框。

        Args:
            body: 按执行顺序排列的指令序列。
            controls: 从外层 ``Control`` 节点继承的控制点；省略为无控制。
        """
        for node in body:
            if isinstance(node, Primitive):
                self.emit_primitive(node, controls)
            elif isinstance(node, Load):
                self.emit_load(node, controls)
            elif isinstance(node, Store):
                self.emit_store(node, controls)
            elif isinstance(node, Call):
                self.emit_call(node, controls)
            elif isinstance(node, Repeat):
                count = node.count
                if count > 1024 and count & (count - 1) == 0:
                    label = f"\\times 2^{{{count.bit_length() - 1}}}"
                else:
                    label = f"\\times {count}"
                self.emit_grouped(node, controls, label)
            elif isinstance(node, Adjoint):
                self.emit_grouped(node, controls, "\\dagger")
            elif isinstance(node, Control):
                bits = tuple(
                    (row, (node.value >> offset) & 1)
                    for offset, row in enumerate(self.rows(node.register))
                )
                self.emit(node.body, controls + bits)
            else:
                raise ValidationError(f"quantikz 暂不支持的结构节点：{type(node).__name__}")

    def render(self) -> str:
        """渲染完整的 quantikz 环境，含分组框的 overlay 节点。

        Returns:
            str: 完整的 quantikz LaTeX 源码；有分组时附加 overlay TikZ 节点。
        """
        if self.module.body is None:
            raise ValidationError("开放声明没有本体；请先绑定或直接绘制实现模块")
        self.emit(self.module.body)
        lines: list[str] = []
        for row, (register, bit) in enumerate(self.wires):
            label = (
                f"\\lstick{{$\\mathtt{{{_escape(register)}}}_{{{bit}}}$}}"
                if self.wire_labels
                else ""
            )
            cells = [cell or "\\qw" for cell in (column[row] for column in self.columns)]
            lines.append(" & ".join([label, *cells]) + r" \\")
        prefix = (
            "\\begin{quantikz}[remember picture, every matrix/.append style={name="
            + self.name
            + "}]\n"
        )
        body = prefix + "\n".join(lines) + "\n\\end{quantikz}"
        if not self.groups:
            return body
        nodes: list[str] = []
        for used, label in self.groups:
            anchors = "".join(
                f"({self.name}-{row + 1}-{column + 1})" for column, row in used
            )
            nodes.append(
                "\\node[draw=black!55, dashed, rounded corners=2pt, inner sep=3pt, "
                f"fit={anchors}, "
                f"label={{[font=\\footnotesize\\color{{black!65}}]above:${label}$}}] {{}};"
            )
        return (
            body
            + "\n\\begin{tikzpicture}[remember picture, overlay]\n"
            + "\n".join(nodes)
            + "\n\\end{tikzpicture}"
        )


def quantikz(
    source: Program | Operation,
    *,
    wire_labels: bool = True,
    name: str | None = None,
) -> str:
    """把 Operation 或 Program 的入口模块绘制成 quantikz 代码。

    Call 绘制成带模块名的盒子（开放 oracle 槽位用虚线），Repeat/Adjoint
    的组框通过 remember picture overlay 标注次数与 dagger，Control 展开为
    控制点；QRAM Load 绘制成资源盒子。结构保持与 RIR 一致，不做任何展开
    或内联。name 固定矩阵节点名，便于同一文档中放置多幅线路。

    Args:
        source: ``Program``、``Operation`` 或带 ``operation`` 属性的输入模型对象。
        wire_labels: 是否在每条量子线左侧绘制寄存器名与位下标标签。
        name: TikZ 矩阵节点名，供同一文档放置多幅线路；省略时按模块名自动生成。

    Returns:
        str: 入口模块的 quantikz LaTeX 源码。
    """
    if isinstance(source, Program):
        program = source
    elif hasattr(source, "program"):
        program = source.program()
    elif hasattr(source, "operation"):
        program = source.operation.program()
    else:
        raise ValidationError(f"quantikz 无法识别的输入类型：{type(source).__name__}")
    if name is None:
        safe = "".join(ch for ch in program.main.name if ch.isalnum()) or "circuit"
        name = f"pq{safe}{next(_name_counter)}"
    return _Renderer(program, wire_labels, name).render()
