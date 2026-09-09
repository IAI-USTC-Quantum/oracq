"宿主读出计划；量子 RIR 保持可组合，末端测量和重置单独描述。"

from dataclasses import dataclass, replace

from pyqecclang.infrastructure.backends import export_originir
from pyqecclang.infrastructure.ir import ValidationError


@dataclass(frozen=True)
class ReadoutAction:
    kind: str
    register: str


def export_with_readout(program, actions=()):
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
