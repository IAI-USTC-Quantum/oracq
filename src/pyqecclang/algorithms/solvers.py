"""可注入 QLSS/HamSim 的 QODE/QPDE 生成器，侧重组装范式。"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..builder import Builder
from ..combinators import lcu, projector, tensor, truncated_shift
from ..ir import Bits, ValidationError, fuse
from ..library import BlockEncoding, _name, identity
from ..oracles import StateOracle, StatePreparation, annotate, invoke, resources_for


def extend_initial(prep, extra_width):
    b = Builder(
        _name("extend_initial", prep.operation, extra_width),
        {"target": Bits(prep.width + extra_width), "work": Bits(prep.work_width)},
        resources_for(("initial", prep.operation)),
    )
    invoke(b, prep.operation, "initial", target=b["target"][: prep.width], work=b["work"])
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def select_subspace(state: StateOracle, output_width, high_value=0, *, label="selection"):
    extra = state.width - output_width
    if extra < 0 or not 0 <= high_value < 1 << extra:
        raise ValidationError("输出子空间布局无效")
    b = Builder(
        _name(label, state.operation, output_width, high_value),
        {"target": Bits(output_width), "signal": Bits(state.signal_qubits + extra)},
        resources_for(("state", state.operation)),
        attributes={
            "algorithm": label,
            "selected_high_value": high_value,
            "validation_stage": "paradigm",
            "success_condition": "signal == 0",
        },
    )
    old_signal, high = b["signal"][: state.signal_qubits], b["signal"][state.signal_qubits :]
    invoke(b, state.operation, "state", target=fuse(b["target"], high), signal=old_signal)
    for bit in range(extra):
        if (high_value >> bit) & 1:
            b.x(high[bit])
    return StateOracle(b.finish())


def make_euler_history_qode(qlss, *, steps=2):
    if type(steps) is not int or steps < 1:
        raise ValidationError("时间步数必须为正整数")

    def generate(generator: BlockEncoding, initial: StatePreparation, final_time):
        if generator.width != initial.width or final_time <= 0:
            raise ValidationError("QODE 输入布局或时间无效")
        nt = steps.bit_length()
        dt = final_time / steps
        q = projector(nt, range(1, steps + 1))
        shift = truncated_shift(nt, steps)
        c = lcu(
            [
                (1, identity(nt + generator.width)),
                (-dt, tensor(q, generator)),
                (-1, tensor(shift, identity(generator.width))),
            ]
        )
        rhs = extend_initial(initial, nt)
        history = qlss(c, rhs)
        return select_subspace(history, generator.width, steps, label="qode_final_time")

    return generate


@dataclass(frozen=True)
class DiscretePDE:
    generator: BlockEncoding
    initial: StatePreparation
    label: str = "linear_pde"


def make_qpde(qode, discretizer=lambda problem: problem):
    def generate(problem, final_time):
        discrete = discretizer(problem)
        result = qode(discrete.generator, discrete.initial, final_time)
        return result

    return generate


def apply_be_to_state(a, prep):
    if a.width != prep.width:
        raise ValidationError("算子与制备的目标宽度不符")
    b = Builder(
        _name("apply_be_state", a.operation, prep.operation),
        {"target": Bits(a.width), "signal": Bits(a.signal_qubits + prep.work_width)},
        resources_for(("a", a.operation), ("prep", prep.operation)),
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["signal"][a.signal_qubits :])
    invoke(b, a.operation, "a", target=b["target"], signal=b["signal"][: a.signal_qubits])
    return StateOracle(b.finish())


def make_lchs_qode(ham_sim, times, weights):
    times, weights = tuple(times), tuple(weights)
    if len(times) != len(weights):
        raise ValidationError("LCHS 离散时间与权重长度不同")

    def generate(generator, initial, final_time):
        terms = []
        for weight, time in zip(weights, times, strict=True):
            op = ham_sim(generator, time * final_time)
            terms.append((weight, BlockEncoding(annotate(op, "block_encoding", be_alpha=1.0))))
        return apply_be_to_state(lcu(terms), initial)

    return generate


def trotter_hamsim(terms, final_time, *, steps=2):
    terms = tuple((float(c), word) for c, word in terms)
    if not terms or steps < 1:
        raise ValidationError("Trotter 需要非空项和正步数")
    width = len(terms[0][1])
    if any(len(word) != width for _, word in terms):
        raise ValidationError("Pauli 项宽度不同")
    b = Builder(
        _name("trotter", terms, final_time, steps),
        {"target": Bits(width), "signal": Bits(0)},
        attributes={"algorithm": "trotter_hamsim", "validation_stage": "paradigm"},
    )
    with b.repeat(steps):
        for coefficient, word in terms:
            active = [i for i, char in enumerate(word) if char != "I"]
            if not active:
                b.global_phase(-coefficient * final_time / steps)
                continue
            for bit in active:
                if word[bit] == "Y":
                    b.gate("phase", b["target"][bit], -math.pi / 2)
                if word[bit] in "XY":
                    b.h(b["target"][bit])
            for bit in active[:-1]:
                b.xor(b["target"][bit], b["target"][active[-1]])
            b.rz(b["target"][active[-1]], 2 * coefficient * final_time / steps)
            for bit in reversed(active[:-1]):
                b.xor(b["target"][bit], b["target"][active[-1]])
            for bit in reversed(active):
                if word[bit] in "XY":
                    b.h(b["target"][bit])
                if word[bit] == "Y":
                    b.gate("phase", b["target"][bit], math.pi / 2)
    return b.finish()


def make_schrodingerisation_qode(embedding, ham_sim, *, extra_width=1):
    """显式保留 Hamiltonian lift 的实现边界，随后组装演化和物理通道。"""

    def generate(generator, initial, final_time):
        hamiltonian = embedding(generator)
        if hamiltonian.width != generator.width + extra_width:
            raise ValidationError("Schrodingerisation lift 的寄存器宽度不符")
        evolution = ham_sim(hamiltonian, final_time)
        encoded = BlockEncoding(annotate(evolution, "block_encoding", be_alpha=1.0))
        lifted_state = apply_be_to_state(encoded, extend_initial(initial, extra_width))
        return select_subspace(
            lifted_state, generator.width, 0, label="schrodingerisation_physical_channel"
        )

    return generate
