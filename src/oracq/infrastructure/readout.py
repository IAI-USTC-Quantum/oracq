"宿主读出计划；量子 RIR 保持可组合，末端测量和重置单独描述。"

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

from oracq.infrastructure.backends import OriginIRArtifact, export_originir
from oracq.infrastructure.ir import Program, ValidationError


@dataclass(frozen=True)
class ReadoutAction:
    """一条末端宿主读出动作的描述。

    Attributes:
        kind: 动作类别，``measure`` 或 ``reset``。
        register: 作用的公开寄存器名。
    """

    kind: str
    register: str


def export_with_readout(
    program: Program, actions: Sequence[ReadoutAction] = ()
) -> OriginIRArtifact:
    """导出闭合程序的 OriginIR，并按序在末端追加宿主读出动作。

    Args:
        program: 已闭合、待导出的程序。
        actions: 依次执行的 ``ReadoutAction`` 序列；为空时直接返回导出结果。

    Returns:
        OriginIRArtifact: 追加动作后的导出产物；无动作时为原样导出结果。

    Raises:
        ValidationError: 动作类别不是测量或重置，或寄存器名不在导出产物中。

    仅在提供动作时才延迟导入 ``uniqc``，把模块化文本展平为扩展 OriginIR
    （宿主动态解析器不接收 DEF）；测量结果写入顺序编号的经典位，
    并相应改写 CREG 行。
    """
    artifact = export_originir(program)
    if not actions:
        return artifact
    # 当前 UnifiedQuantum 动态解析器不接收 DEF，只有显式宿主读出适配才展平。
    from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

    parser = OriginIR_BaseParser()
    parser.parse(artifact.text)
    artifact = replace(artifact, text=parser.to_extended_originir())
    tail, count = [], 0
    for action in actions:
        if action.kind not in {"measure", "reset"} or action.register not in artifact.registers:
            raise ValidationError("无效的末端读出动作")
        for qubit in artifact.registers[action.register]:
            if action.kind == "measure":
                tail.append(f"MEASURE q[{qubit}], c[{count}]")
                count += 1
            else:
                tail.append(f"RESET q[{qubit}]")
    lines = artifact.text.splitlines()
    lines[lines.index("CREG 0")] = f"CREG {count}"
    return replace(artifact, text="\n".join([*lines, *tail]) + "\n")
