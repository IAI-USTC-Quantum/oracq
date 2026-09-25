"""quantikz circuit export of RIR modules: module structure is preserved, without expanding calls and repeats."""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence
from typing import cast

from oracq.infrastructure.builder import Operation
from oracq.infrastructure.ir import (
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
    """Escape underscores in names for LaTeX math mode."""
    return name.replace("_", r"\_")


def _angle_text(angle: float) -> str:
    """Render small integer multiples of pi as fraction-bar form, the rest as decimals."""
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
    """Split a list of row numbers into consecutive runs, preserving the given order."""
    result: list[list[int]] = []
    for row in rows:
        if result and row == result[-1][-1] + 1:
            result[-1].append(row)
        else:
            result.append([row])
    return result


class _Renderer:
    """Lay the entry module out node by node into quantikz matrix columns, and register grouping frames."""

    def __init__(self, program: Program, wire_labels: bool, name: str) -> None:
        """Bind the program, wire-label option and module name, and initialize the quantum wire table."""
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
        """Expand a view into the list of matrix rows it covers.

        Args:
            ref: The register view to expand.

        Returns:
            list[int]: The quantum-wire row numbers covered by the view, ordered from least to most significant bit within each span.
        """
        result: list[int] = []
        for span in ref.parts:
            result.extend(
                self.row_of[(span.register, span.start + offset)] for offset in range(span.width)
            )
        return result

    def new_column(self) -> list[str]:
        """Append one column of empty cells and return that column.

        Returns:
            list[str]: The newly appended column; each element is the cell of one quantum wire.
        """
        self.columns.append([""] * len(self.wires))
        return self.columns[-1]

    def place_controls(
        self, column: list[str], target_rows: Sequence[int], controls: tuple[tuple[int, int], ...]
    ) -> None:
        """Place control points in the column, pointing at the first row of each target segment.

        Args:
            column: The matrix column to write control-point marks into.
            target_rows: The set of target row numbers the controls point at; the minimum row is the pointed end.
            controls: Sequence of (row number, active bit value) tuples; a bit value of 0 draws a hollow anti-control point.
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
        """Draw each consecutive run of row numbers as a box gate, adding control points as needed.

        Args:
            rows: The row numbers covered by the box; each consecutive run becomes its own box occupying one column.
            label: The LaTeX label inside the box.
            controls: Control points appended to each box.
            dashed: When true the box uses a dashed style, marking unbound structures such as open declarations.
        """
        for run in _runs(rows):
            column = self.new_column()
            style = ",style={dashed}" if dashed else ""
            column[run[0]] = f"\\gate[wires={len(run)}{style}]{{{label}}}"
            if controls:
                self.place_controls(column, run, controls)

    def emit_primitive(self, node: Primitive, controls: tuple[tuple[int, int], ...]) -> None:
        """Emit one primitive gate: single-qubit gates, CNOT, SWAP, add-constant or global phase.

        Args:
            node: The primitive instruction node to draw.
            controls: Control points inherited from outer ``Control`` nodes.
        """
        op = node.op
        if op in _GATE_LABELS or op in _ROTATION_LABELS:
            label = _GATE_LABELS.get(op)
            if op in _ROTATION_LABELS:
                # Validation guarantees the rotation gate angle exists.
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
                # Validation guarantees the gphase angle exists.
                f"e^{{i{_angle_text(cast(float, node.angle))}}}",
                controls,
            )
            return
        raise ValidationError(f"quantikz does not yet support the gate operation: {op}")

    def emit_load(self, node: Load, controls: tuple[tuple[int, int], ...]) -> None:
        """Draw a QRAM query as a resource box.

        Args:
            node: The QRAM read instruction to draw.
            controls: Control points inherited from outer ``Control`` nodes.
        """
        rows = self.rows(node.address) + self.rows(node.data)
        self.box(rows, f"\\mathrm{{QRAM}}_{{\\mathtt{{{_escape(node.resource)}}}}}", controls)

    def emit_store(self, node: Store, controls: tuple[tuple[int, int], ...]) -> None:
        """Draw a QRAM random write as a resource box with superscript ``w``.

        Args:
            node: The QRAM write instruction to draw.
            controls: Control points inherited from outer ``Control`` nodes.
        """
        rows = self.rows(node.address) + self.rows(node.data)
        self.box(
            rows,
            f"\\mathrm{{QRAM}}^{{\\mathrm{{w}}}}_{{\\mathtt{{{_escape(node.resource)}}}}}",
            controls,
        )

    def emit_call(self, node: Call, controls: tuple[tuple[int, int], ...]) -> None:
        """Draw a module call as a box labeled with the module name; open declarations use dashes.

        Args:
            node: The module call instruction to draw.
            controls: Control points inherited from outer ``Control`` nodes.
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
        """Emit the instructions inside the group, then register one annotated grouping frame over the cells they occupy.

        Args:
            node: The ``Repeat`` or ``Adjoint`` structure node whose body is emitted.
            controls: Control points inherited from outer ``Control`` nodes.
            label: The LaTeX annotation above the grouping frame, such as the repeat count or dagger.
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
            # A single box collapses into a superscript mark (the DB^\dagger form is more compact).
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
        """Dispatch emission by node type; structural nodes recurse or collapse into grouping frames.

        Args:
            body: The instruction sequence in execution order.
            controls: Control points inherited from outer ``Control`` nodes; no controls when omitted.
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
                raise ValidationError(f"quantikz does not yet support the structural node: {type(node).__name__}")

    def render(self) -> str:
        """Render the complete quantikz environment, including the overlay nodes of grouping frames.

        Returns:
            str: The complete quantikz LaTeX source; overlay TikZ nodes are appended when groupings exist.
        """
        if self.module.body is None:
            raise ValidationError("open declarations have no body; bind one first, or draw an implementation module directly")
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
    """Draw the entry module of an Operation or Program as quantikz code.

    Calls are drawn as boxes carrying the module name (open oracle slots use
    dashes), Repeat/Adjoint group frames annotate the count and dagger via
    remember picture overlays, and Controls expand into control points;
    QRAM Loads are drawn as resource boxes. The structure stays identical
    to RIR, with no expansion or inlining. name fixes the matrix node name
    so several circuits can be placed in the same document.

    Args:
        source: A ``Program``, an ``Operation``, or an input-model object with an ``operation`` attribute.
        wire_labels: Whether to draw register name and bit index labels on the left of each quantum wire.
        name: TikZ matrix node name, for placing several circuits in one document; generated from the module name when omitted.

    Returns:
        str: The quantikz LaTeX source of the entry module.
    """
    if isinstance(source, Program):
        program = source
    elif hasattr(source, "program"):
        program = source.program()
    elif hasattr(source, "operation"):
        program = source.operation.program()
    else:
        raise ValidationError(f"quantikz cannot recognize the input type: {type(source).__name__}")
    if name is None:
        safe = "".join(ch for ch in program.main.name if ch.isalnum()) or "circuit"
        name = f"pq{safe}{next(_name_counter)}"
    return _Renderer(program, wire_labels, name).render()
