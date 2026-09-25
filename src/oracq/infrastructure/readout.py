"Host readout plan; the quantum RIR stays composable, with terminal measurement and reset described separately."

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

from oracq.infrastructure.backends import OriginIRArtifact, export_originir
from oracq.infrastructure.ir import Program, ValidationError


@dataclass(frozen=True)
class ReadoutAction:
    """Description of one terminal host readout action.

    Attributes:
        kind: action category, ``measure`` or ``reset``.
        register: public register name the action applies to.
    """

    kind: str
    register: str


def export_with_readout(
    program: Program, actions: Sequence[ReadoutAction] = ()
) -> OriginIRArtifact:
    """Export a closed program to OriginIR and append host readout actions in order at the end.

    Args:
        program: closed program to export.
        actions: sequence of ``ReadoutAction`` steps executed in order; when empty the export result is returned directly.

    Returns:
        OriginIRArtifact: export artifact with the actions appended; the plain export result when there are no actions.

    Raises:
        ValidationError: an action category is neither measurement nor reset, or a register name is absent from the export artifact.

    ``uniqc`` is imported lazily only when actions are provided, flattening the
    modular text into extended OriginIR (the host dynamic parser does not
    accept DEF); measurement results go into sequentially numbered classical
    bits, and the CREG line is rewritten accordingly.
    """
    artifact = export_originir(program)
    if not actions:
        return artifact
    # The current UnifiedQuantum dynamic parser does not accept DEF; flattening happens only for explicit host readout adaptation.
    from uniqc.compile.originir.originir_base_parser import OriginIR_BaseParser

    parser = OriginIR_BaseParser()
    parser.parse(artifact.text)
    artifact = replace(artifact, text=parser.to_extended_originir())
    tail, count = [], 0
    for action in actions:
        if action.kind not in {"measure", "reset"} or action.register not in artifact.registers:
            raise ValidationError("invalid terminal readout action")
        for qubit in artifact.registers[action.register]:
            if action.kind == "measure":
                tail.append(f"MEASURE q[{qubit}], c[{count}]")
                count += 1
            else:
                tail.append(f"RESET q[{qubit}]")
    lines = artifact.text.splitlines()
    lines[lines.index("CREG 0")] = f"CREG {count}"
    return replace(artifact, text="\n".join([*lines, *tail]) + "\n")
