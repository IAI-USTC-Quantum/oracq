"""Toffoli+Clifford+T+QRAM strict netlist export.

Runs one more U3 classification pass over the output of
``basis.lower_toffoli_u3_cz``: rotations whose angles are integer multiples
of π/4 fall into exact Clifford+T atoms (H/S/SDG/T/TDG/Z/X/Y), while the
rest stay as ``RZ``/``RY`` rotation atoms pending synthesis; module calls,
Repeat symbolic structures and QRAM lines remain unexpanded. Atom set:
TOFFOLI, CZ, H, S, SDG, T, TDG, X, Y, Z, RZ, RY, QRAM. Angle classification
shares the same implementation with the ``estimate`` module, keeping netlist
counts and resource estimation consistent item by item (the tests
cross-check them line by line)."""

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
    "Strict netlist: text + per-atom counts + layout information."

    text: str
    counts: Counter[str]
    registers: dict[str, tuple[int, ...]]
    resources: dict[str, str]
    workspace_qubits: tuple[int, ...]
    qram_queries: Counter[str]
    qram_writes: Counter[str] = field(default_factory=Counter)


def _emit_rz(lines: list[str], counts: Counter[str], target: str, angle: float) -> None:
    """Emit the atom sequence of an RZ angle to the netlist per the ``classify_rz`` classification."""
    atoms, rotations = classify_rz(angle)
    for atom, n in sorted(atoms.items()):
        lines.extend(f"{_ATOM_NAMES[atom]} {target}" for _ in range(n))
    for axis, value in rotations:
        lines.append(f"{axis.upper()} {target}, ({value!r})")
    counts.update(atoms)
    counts.update(axis for axis, _ in rotations)


def _emit_ry(lines: list[str], counts: Counter[str], target: str, angle: float) -> None:
    """Emit the atom sequence of an RY angle to the netlist per the ``classify_ry`` classification."""
    atoms, rotations = classify_ry(angle)
    if rotations:
        lines.append(f"RY {target}, ({rotations[0][1]!r})")
        counts["ry"] += 1
        return
    if not atoms:
        return  # Identity: RY(0) emits nothing
    # RY(θ) = S·H·RZ(θ)·H·S†: S/H sandwich + exact rz sequence
    lines.append(f"S {target}")
    lines.append(f"H {target}")
    counts.update({"s": 1, "h": 2, "sdg": 1})
    _emit_rz(lines, counts, target, angle)
    lines.append(f"H {target}")
    lines.append(f"SDG {target}")


def _emit_u3(
    lines: list[str], counts: Counter[str], target: str, theta: float, phi: float, lam: float
) -> None:
    """Emit the atom sequence of U3 angles to the netlist per the ``classify_u3`` classification."""
    atoms, rotations = classify_u3(theta, phi, lam)
    if not rotations and set(atoms) <= {"h", "x", "y"} and sum(atoms.values()) == 1:
        atom = next(iter(atoms))
        lines.append(f"{_ATOM_NAMES[atom]} {target}")
        counts[atom] += 1
        return
    # Same Euler order as classify_u3: RZ(λ) → RY(θ) → RZ(φ)
    _emit_rz(lines, counts, target, lam)
    _emit_ry(lines, counts, target, theta)
    _emit_rz(lines, counts, target, phi)


def lower_strict(artifact: OriginIRArtifact) -> StrictArtifact:
    """Further lower a Toffoli+U3+CZ netlist into a strict atom netlist.

    ``U3`` lines are rewritten by the angle classification shared with
    ``estimate``: angles that are integer multiples of π/4 emit exact
    Clifford+T atoms, while the rest stay as ``RZ``/``RY`` rotation atoms
    pending synthesis; every other line is kept verbatim, while atom counts,
    QRAM queries (``ram_`` lines) and random writes (``QRAMWRITE`` lines)
    are accumulated line by line.

    Args:
        artifact: The netlist produced by ``lower_toffoli_u3_cz``.

    Returns:
        StrictArtifact: The strict netlist text, per-atom counts and original layout information.

    Raises:
        ValidationError: A ``U3`` line does not carry three angle parameters.
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
                raise ValidationError("strict lowering encountered an invalid U3 line: " + line)
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
    """Compile to a module-preserving netlist at the Toffoli+Clifford+T(+rotations pending synthesis)+QRAM level.

    Args:
        program: The closed RIR program to export.

    Returns:
        StrictArtifact: The netlist text, together with Toffoli/CZ gate counts, QRAM query
        and write statistics, and the register and resource mappings.
    """
    return lower_strict(lower_toffoli_u3_cz(export_originir(program)))
