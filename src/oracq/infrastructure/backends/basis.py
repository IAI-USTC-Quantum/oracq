"在 OriginIR DEF 边界内降低门集；模块调用和 QRAM 声明保持不展开。"

from __future__ import annotations

import math
import re

from oracq.infrastructure.backends.originir import OriginIRArtifact, export_originir
from oracq.infrastructure.ir import Program, ValidationError

_PI = math.pi


def export_toffoli_u3_cz(program: Program) -> OriginIRArtifact:
    """导出 Toffoli+U3+CZ 门集上的 OriginIR-ext 网表。

    等价于 ``lower_toffoli_u3_cz(export_originir(program))``；模块调用和
    QRAM 声明保持不展开。

    Args:
        program: 待编译的封闭 RIR 程序。

    Returns:
        OriginIRArtifact: 降低后的网表，``workspace_qubits`` 含新增的
        ``pb_work`` 辅助比特编号。
    """
    return lower_toffoli_u3_cz(export_originir(program))


def lower_toffoli_u3_cz(artifact: OriginIRArtifact) -> OriginIRArtifact:
    """把 OriginIR-ext 网表中的门降低到 Toffoli、U3 与 CZ。

    X 型多控门（含 CNOT，SWAP 展开为三次交换）按 ``mcx`` 配方用 Toffoli
    梯子实现，受控单比特门按 ``controlled_u3`` 配方实现，两者共享按最大
    控制数分配的 ``pb_work`` 辅助比特池；池位追加进 ``QINIT``、``DEF``
    形参与 ``m_`` 调用行。``DAGGER``、QRAM 与 ``CREG`` 行原样保留。

    Args:
        artifact: ``export_originir`` 产出的 OriginIR-ext 网表。

    Returns:
        OriginIRArtifact: 降低后的网表；``workspace_qubits`` 追加新增的
        辅助比特编号，``registers`` 与 ``resources`` 原样透传。

    Raises:
        ValidationError: 遇到参数表不支持的目标门。
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
        """发射单比特 ``U3`` 门文本行。"""
        return f"U3 {t}, ({float(theta)!r}, {float(phi)!r}, {float(lam)!r})"

    def h(t: str) -> str:
        """发射 ``H`` 门的 ``U3`` 等价文本行。"""
        return u3(t, _PI / 2, 0, _PI)

    def phase(t: str, angle: float) -> str:
        """发射相位旋转门的 ``U3`` 等价文本行。"""
        return u3(t, 0, 0, angle)

    def cx(c: str, t: str) -> list[str]:
        """用 H 环绕 CZ 实现 ``CNOT``。"""
        return [h(t), f"CZ {c}, {t}", h(t)]

    def mcx(controls: list[str], target: str) -> list[str]:
        """用 Toffoli 梯子（含 ``pb_work`` 辅助比特池）实现多控 X。"""
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
        """用 Toffoli 梯子加相位分解实现受控 ``U3``。"""
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
                raise ValidationError("不支持的目标门集 lowering 输入：" + gate)
            result += controlled_u3(controls, bits[0], *parameters[gate])
    return OriginIRArtifact(
        "\n".join(result) + "\n",
        artifact.registers,
        artifact.resources,
        artifact.workspace_qubits + tuple(range(total, total + count)),
    )
