"Lower the gate set within OriginIR DEF boundaries; module calls and QRAM declarations remain unexpanded."

from __future__ import annotations

import math
import re

from oracq.infrastructure.backends.originir import OriginIRArtifact, export_originir
from oracq.infrastructure.ir import Program, ValidationError

_PI = math.pi


def export_toffoli_u3_cz(program: Program) -> OriginIRArtifact:
    """Export an OriginIR-ext netlist over the Toffoli+U3+CZ gate set.

    Equivalent to ``lower_toffoli_u3_cz(export_originir(program))``; module
    calls and QRAM declarations remain unexpanded.

    Args:
        program: The closed RIR program to compile.

    Returns:
        OriginIRArtifact: The lowered netlist; ``workspace_qubits`` includes the newly added
        ``pb_work`` ancilla bit indices.
    """
    return lower_toffoli_u3_cz(export_originir(program))


def lower_toffoli_u3_cz(artifact: OriginIRArtifact) -> OriginIRArtifact:
    """Lower the gates of an OriginIR-ext netlist to Toffoli, U3 and CZ.

    Multi-controlled X-type gates (CNOT included, SWAP expanded into three
    swaps) are implemented with a Toffoli ladder per the ``mcx`` recipe,
    and controlled single-qubit gates per the ``controlled_u3`` recipe;
    both share the ``pb_work`` ancilla pool allocated by the maximum
    control count; pool bits are appended to ``QINIT``, ``DEF`` formals and
    ``m_`` call lines. ``DAGGER``, QRAM and ``CREG`` lines are kept
    verbatim.

    Args:
        artifact: The OriginIR-ext netlist produced by ``export_originir``.

    Returns:
        OriginIRArtifact: The lowered netlist; ``workspace_qubits`` appends the newly added
        ancilla bit indices, while ``registers`` and ``resources`` pass through unchanged.

    Raises:
        ValidationError: A target gate unsupported by the parameter table was encountered.
    """
    lines = artifact.text.splitlines()
    maximum = 0
    for line in lines:
        head, _, tail = line.partition(" controlled_by (")
        controls = tail[:-1].split(", ") if tail else []
        gate = head.split(" ", 1)[0]
        if gate in {"CNOT", "CZ"}:
            controls.append("implicit")
        if gate not in {"DEF", "ENDDEF", "QRAMWRITE"} and not gate.startswith(("m_", "ram_")):
            maximum = max(maximum, len(controls))
    count = max(0, maximum - 1)
    total = int(next(line.split()[1] for line in lines if line.startswith("QINIT ")))
    pool = [f"pb_work[{i}]" for i in range(count)]
    result, inside = [], False

    def u3(t: str, theta: float, phi: float, lam: float) -> str:
        """Emit one single-qubit ``U3`` gate text line."""
        return f"U3 {t}, ({float(theta)!r}, {float(phi)!r}, {float(lam)!r})"

    def h(t: str) -> str:
        """Emit the ``U3``-equivalent text line of the ``H`` gate."""
        return u3(t, _PI / 2, 0, _PI)

    def phase(t: str, angle: float) -> str:
        """Emit the ``U3``-equivalent text line of a phase rotation gate."""
        return u3(t, 0, 0, angle)

    def cx(c: str, t: str) -> list[str]:
        """Implement ``CNOT`` with an H-sandwiched CZ."""
        return [h(t), f"CZ {c}, {t}", h(t)]

    def mcx(controls: list[str], target: str) -> list[str]:
        """Implement multi-controlled X with a Toffoli ladder, using the ``pb_work`` ancilla pool."""
        n = len(controls)
        if n == 0:
            return [u3(target, _PI, 0, _PI)]
        if n == 1:
            return cx(controls[0], target)
        if n == 2:
            return [f"TOFFOLI {controls[0]}, {controls[1]}, {target}"]
        ladder = [f"TOFFOLI {controls[0]}, {controls[1]}, {pool[0]}"]
        ladder += [f"TOFFOLI {pool[i - 2]}, {controls[i]}, {pool[i - 1]}" for i in range(2, n - 1)]
        return (
            ladder + [f"TOFFOLI {pool[n - 3]}, {controls[-1]}, {target}"] + list(reversed(ladder))
        )

    def controlled_u3(
        controls: list[str],
        target: str,
        theta: float,
        phi: float,
        lam: float,
        global_angle: float = 0,
    ) -> list[str]:
        """Implement a controlled ``U3`` with a Toffoli ladder plus phase decomposition."""
        if not controls:
            output = [u3(target, theta, phi, lam)]
            if global_angle:
                output += [
                    phase(target, global_angle),
                    *mcx([], target),
                    phase(target, global_angle),
                    *mcx([], target),
                ]
            return output
        ladder: list[str] = []
        if len(controls) == 1:
            c = controls[0]
        else:
            ladder = [f"TOFFOLI {controls[0]}, {controls[1]}, {pool[0]}"]
            ladder += [
                f"TOFFOLI {pool[i - 2]}, {controls[i]}, {pool[i - 1]}"
                for i in range(2, len(controls))
            ]
            c = pool[len(controls) - 2]
        middle = [phase(c, (lam + phi) / 2 + global_angle), phase(target, (lam - phi) / 2)]
        middle += cx(c, target)
        middle += [u3(target, -theta / 2, 0, -(phi + lam) / 2)]
        middle += cx(c, target)
        middle += [u3(target, theta / 2, phi, 0)]
        return ladder + middle + list(reversed(ladder))

    for line in lines:
        if line.startswith("QINIT "):
            result.append(f"QINIT {total + count}")
            continue
        if line.startswith("DEF "):
            inside = True
            if count:
                line = line[:-1] + (", " if not line.endswith("()") else "") + f"pb_work[{count}])"
            result.append(line)
            continue
        if line == "ENDDEF":
            inside = False
            result.append(line)
            continue
        if line.startswith("m_"):
            if count:
                extra = pool if inside else [f"q[{i}]" for i in range(total, total + count)]
                line = (
                    line[:-1] + (", " if not line.endswith("()") else "") + ", ".join(extra) + ")"
                )
            result.append(line)
            continue
        if line in {"DAGGER", "ENDDAGGER"} or line.startswith(("QRAMDECL ", "QRAMWRITE ", "CREG ", "ram_")):
            result.append(line)
            continue
        head, _, tail = line.partition(" controlled_by (")
        controls = tail[:-1].split(", ") if tail else []
        gate, _, operands = head.partition(" ")
        bits = re.findall(r"[A-Za-z_][A-Za-z0-9_]*\[\d+\]", operands)
        angle_match = re.search(r"\(([^()]*)\)", operands)
        angle = float(angle_match.group(1)) if angle_match else 0.0
        if gate == "CNOT":
            result += mcx(controls + [bits[0]], bits[1])
        elif gate == "SWAP":
            result += mcx(controls + [bits[0]], bits[1])
            result += mcx(controls + [bits[1]], bits[0])
            result += mcx(controls + [bits[0]], bits[1])
        elif gate == "X":
            result += mcx(controls, bits[0])
        elif gate == "CZ":
            result += controlled_u3(controls + [bits[0]], bits[1], 0, 0, _PI)
        else:
            parameters = {
                "H": (_PI / 2, 0, _PI, 0),
                "Y": (_PI, _PI / 2, _PI / 2, 0),
                "Z": (0, 0, _PI, 0),
                "S": (0, 0, _PI / 2, 0),
                "T": (0, 0, _PI / 4, 0),
                "U1": (0, 0, angle, 0),
                "RY": (angle, 0, 0, 0),
                "RX": (angle, -_PI / 2, _PI / 2, 0),
                "RZ": (0, 0, angle, -angle / 2),
            }
            if gate not in parameters:
                raise ValidationError("unsupported target gate set lowering input: " + gate)
            result += controlled_u3(controls, bits[0], *parameters[gate])
    return OriginIRArtifact(
        "\n".join(result) + "\n",
        artifact.registers,
        artifact.resources,
        artifact.workspace_qubits + tuple(range(total, total + count)),
    )
