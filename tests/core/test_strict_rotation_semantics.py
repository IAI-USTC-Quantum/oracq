"""Independent matrix replay checks strict rotation signs and controlled phases."""

import cmath
import math
import re
import unittest

from oracq import Bits, Builder, export_strict, simulate


def replay(text, n, column):
    state = [complex(i == column) for i in range(1 << n)]
    instructions = []
    stack = []
    for line in text.splitlines():
        gate = line.split()[0]
        if gate == "DAGGER":
            stack.append(instructions)
            instructions = []
        elif gate == "ENDDAGGER":
            instructions = stack.pop() + [(g, q, a, not inv) for g, q, a, inv in reversed(instructions)]
        elif gate in {"H", "S", "SDG", "T", "TDG", "X", "Y", "Z", "CZ", "RZ", "RY"}:
            qubits = tuple(map(int, re.findall(r"v_q\[(\d+)\]", line)))
            angle = float(re.search(r"\(([^)]+)\)", line).group(1)) if gate in {"RZ", "RY"} else None
            instructions.append((gate, qubits, angle, False))
        elif gate not in {"QINIT", "CREG", "DEF", "ENDDEF"} and not gate.startswith("m_"):
            raise AssertionError("unsupported strict instruction: " + line)
    for gate, qubits, angle, inverse in instructions:
        if gate == "CZ":
            state = [-a if all(i & (1 << q) for q in qubits) else a for i, a in enumerate(state)]
            continue
        matrices = {
            "H": ((1 / math.sqrt(2), 1 / math.sqrt(2)), (1 / math.sqrt(2), -1 / math.sqrt(2))),
            "X": ((0, 1), (1, 0)), "Y": ((0, -1j), (1j, 0)), "Z": ((1, 0), (0, -1)),
            "S": ((1, 0), (0, 1j)), "SDG": ((1, 0), (0, -1j)),
            "T": ((1, 0), (0, cmath.exp(1j * math.pi / 4))),
            "TDG": ((1, 0), (0, cmath.exp(-1j * math.pi / 4))),
        }
        if gate == "RZ":
            matrix = ((cmath.exp(-1j * angle / 2), 0), (0, cmath.exp(1j * angle / 2)))
        elif gate == "RY":
            c, s = math.cos(angle / 2), math.sin(angle / 2)
            matrix = ((c, -s), (s, c))
        else:
            matrix = matrices[gate]
        if inverse:
            matrix = tuple(tuple(complex(matrix[j][i]).conjugate() for j in range(2)) for i in range(2))
        bit = 1 << qubits[0]
        updated = state.copy()
        for i in range(1 << n):
            if not i & bit:
                a, b = state[i], state[i | bit]
                updated[i] = matrix[0][0] * a + matrix[0][1] * b
                updated[i | bit] = matrix[1][0] * a + matrix[1][1] * b
        state = updated
    return state


class StrictRotationSemanticsTests(unittest.TestCase):
    def check_matrix(self, program):
        artifact = export_strict(program)
        n = program.main.registers[0].type.width
        phase = None
        for column in range(1 << n):
            actual = replay(artifact.text, n, column)
            expected = simulate(program, initial={"q": column}).statevector()
            for a, e in zip(actual, expected, strict=True):
                if phase is None and abs(e) > 1e-10:
                    phase = a / e
                if phase is not None:
                    self.assertAlmostEqual(a, phase * e, places=12)

    def test_exact_positive_negative_ry_full_matrix(self):
        for k in (-7, -3, -1, 1, 3, 7):
            with self.subTest(k=k):
                b = Builder("ry_exact", {"q": Bits(1)})
                b.ry(b["q"], k * math.pi / 4)
                self.check_matrix(b.finish().program())

    def test_controlled_ry_and_adjoint_preserve_relative_phase(self):
        for k in (-2, -1, 1, 2):
            for inverse in (False, True):
                b = Builder("controlled_ry", {"q": Bits(2)})
                def emit(builder, angle):
                    with builder.control(builder["q"][0], 0):
                        builder.ry(builder["q"][1], angle)
                if inverse:
                    with b.adjoint():
                        emit(b, k * math.pi / 4)
                else:
                    emit(b, k * math.pi / 4)
                with self.subTest(k=k, inverse=inverse):
                    self.check_matrix(b.finish().program())


if __name__ == "__main__":
    unittest.main()
