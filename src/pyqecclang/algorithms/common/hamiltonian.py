"""Hamiltonian 分解、Pauli 演化、Trotter 和有限 Taylor 编码。"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from pyqecclang.algorithms.input_model.block_encoding import lcu, pauli_word
from pyqecclang.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
    requires,
)
from pyqecclang.algorithms.input_model.interfaces import BlockEncodingProtocol, as_block_encoding
from pyqecclang.algorithms.input_model.operators import (
    BlockEncoding,
    _name,
    block_encoding,
    identity,
    product,
)
from pyqecclang.algorithms.input_model.oracles import (
    annotate,
    invoke,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


def trotter_hamsim(terms, final_time, *, steps=2):
    """一阶 Trotter 乘积公式模拟 Pauli 分解的 Hamiltonian 演化。

    对 H = sum_j c_j P_j 按 exp(-i*H*t) ≈ (prod_j exp(-i*c_j*P_j*t/steps))**steps
    合成：每个 Pauli 词经基变换与 CNOT 链折叠为末活跃位上的相位旋转后逆序复原。
    电路主体是单个 Repeat(steps) 块，生成与序列化阶段不按步数展开；
    全 I 的项不消耗量子比特，退化为全局相位。

    Args:
        terms: ``(系数, Pauli 词)`` 序对序列；系数为实数（经 ``float`` 转换），
            Pauli 词为等宽的 I/X/Y/Z 字符串。
        final_time: 总演化时间 t。
        steps: 重复步数，至少为 1；步数越大乘积公式误差越小。

    Returns:
        Operation: 裸酉操作而非块编码；寄存器为 target（词宽度）与零宽 signal，
            整体近似 exp(-i*H*t)。

    Raises:
        ValidationError: 项序列为空、steps 小于 1 或各 Pauli 词宽度不一致。
    """
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


def taylor_hamiltonian(hamiltonian, time, *, degree=2):
    """可闭合的普通 Hamiltonian-function BE；可替换为 QSP/HamSim protocol。"""
    require_instance(hamiltonian, BlockEncoding, "taylor_hamiltonian.H")
    positive_integer(degree, "taylor_hamiltonian.degree", minimum=0)
    finite_real(time, "taylor_hamiltonian.time")
    if degree < 0 or not math.isfinite(time):
        raise ValidationError("Taylor 阶数/演化时间无效")
    powers, current = [(1, identity(hamiltonian.width))], identity(hamiltonian.width)
    for k in range(1, degree + 1):
        current = product(hamiltonian, current)
        powers.append(((-1j * time) ** k / math.factorial(k), current))
    out = lcu(powers)
    return BlockEncoding(
        annotate(
            out.operation,
            "block_encoding",
            be_alpha=out.alpha,
            algorithm="truncated_taylor_hamiltonian_function",
            degree=degree,
            time=float(time),
            correctness="pending",
            success_condition="signal == 0",
        )
    )


@runtime_checkable
class HermitianProtocol(Protocol):
    """宿主声明算符厄米性的访问协议。"""

    @property
    def hermitian(self) -> bool:
        """是否声明为 Hermitian 算符；hamiltonian_simulation 仅接受 True。"""
        ...


@runtime_checkable
class EvolvableProtocol(Protocol):
    """宿主声明算符可给出自身酉演化的访问协议。"""

    def evolution(self, time: float):
        """按演化时长返回该算符的酉 Operation。

        Trotter 路径要求返回的演化无需后选择：除 target 外公开寄存器为零宽。
        """
        ...


@runtime_checkable
class TrotterizableProtocol(Protocol):
    """宿主声明算符可分解为 Trotter 项列表的访问协议。"""

    def trotter_list(self):
        """返回构成 Hamiltonian 的 TrotterTerm 序列；hamiltonian_simulation 要求非空。"""
        ...


@dataclass(frozen=True)
class TrotterTerm:
    """乘积公式的单个 Hamiltonian 项：实系数与可演化算符的组合。

    Attributes:
        coefficient: 项的实系数；Trotter 路径以 ``coefficient * time / steps``
            为时长调用 ``operator.evolution``。
        operator: 满足 EvolvableProtocol 的算符。

    Raises:
        ValidationError: coefficient 非有限实数，或 operator 不满足 EvolvableProtocol。
    """

    coefficient: float
    operator: object

    def __post_init__(self):
        finite_real(self.coefficient, "TrotterTerm.coefficient")
        requires(self.operator, EvolvableProtocol, path="TrotterTerm.operator")


@dataclass(frozen=True)
class PauliOperator:
    """单个 Pauli 词算符；满足 HermitianProtocol 与 EvolvableProtocol。

    Attributes:
        word: 非空的 I/X/Y/Z 字符串，宽度不超过 64。
        hermitian: 恒为 True；即 HermitianProtocol 声明。

    Raises:
        ValidationError: word 不是非空 I/X/Y/Z 字符串，或宽度超过 64。
    """

    word: str
    hermitian = True

    def __post_init__(self):
        if (
            not isinstance(self.word, str)
            or not self.word
            or any(c not in "IXYZ" for c in self.word)
        ):
            raise ValidationError("Pauli word 必须是非空 I/X/Y/Z 字符串")
        positive_integer(len(self.word), "PauliOperator.width", maximum=64)

    def block_encoding(self):
        """返回该 Pauli 词的 alpha=1.0 块编码：逐位单量子比特门加零宽 signal。"""
        return pauli_word(self.word)

    def evolution(self, time):
        """返回 ``exp(-1j*word*time)`` 的精确酉演化（单项单步，无乘积公式误差）。"""
        finite_real(time, "PauliOperator.time")
        return trotter_hamsim(((1.0, self.word),), time, steps=1)


@dataclass(frozen=True)
class PauliHamiltonian:
    """等宽 Pauli 词的实系数线性组合；支持块编码与 Trotter 分解两种访问。

    Attributes:
        terms: ``(系数, Pauli 词)`` 序对的元组；系数为有限实数，所有词等宽。
        hermitian: 恒为 True；即 HermitianProtocol 声明。

    Raises:
        ValidationError: 项列表为空、系数非有限实数、Pauli 词非法或宽度不一致。
    """

    terms: tuple[tuple[float, str], ...]
    hermitian = True

    def __post_init__(self):
        object.__setattr__(self, "terms", tuple(tuple(t) for t in self.terms))
        if not self.terms:
            raise ValidationError("PauliHamiltonian 需要非空项列表")
        for coefficient, word in self.terms:
            finite_real(coefficient, "PauliHamiltonian.coefficient")
            PauliOperator(word)
        if len({len(word) for _, word in self.terms}) != 1:
            raise ValidationError("Pauli 项宽度不一致")

    def block_encoding(self):
        """返回非零系数项的 LCU 块编码；无非零项时退化为零算子块编码。"""
        from pyqecclang.algorithms.input_model.operators import zero

        terms = [(c, pauli_word(w)) for c, w in self.terms if c]
        return lcu(terms) if terms else zero(len(self.terms[0][1]))

    def trotter_list(self):
        """把每个 ``(系数, 词)`` 包装为 TrotterTerm 元组返回。"""
        return tuple(TrotterTerm(c, PauliOperator(w)) for c, w in self.terms)


@dataclass(frozen=True)
class EncodedOperator:
    """矩阵性质的宿主声明；可表示非 Hermitian 算符，不冒称其本身是酉操作。"""

    encoding: BlockEncoding
    hermitian: bool

    def __post_init__(self):
        object.__setattr__(self, "encoding", as_block_encoding(self.encoding))
        if type(self.hermitian) is not bool:
            raise ValidationError("hermitian 需要 bool 声明")

    def block_encoding(self):
        """返回构造时携带的块编码。"""
        return self.encoding


def hamiltonian_simulation(operator, time, *, method="auto", steps=2, qsp=None):
    """当前优先可分解的 Trotter；QSP 需要调用者提供实际实现。"""
    requires(operator, HermitianProtocol, path="HamSim.operator")
    if operator.hermitian is not True:
        raise ValidationError(
            "Hamiltonian simulation 需要 Hermitian；非 Hermitian 动力学请使用 QODE"
        )
    finite_real(time, "HamSim.time")
    if method not in {"auto", "trotter", "qsp"}:
        raise ValidationError("未知 Hamiltonian simulation 方法")
    if method == "auto":
        method = "trotter" if isinstance(operator, TrotterizableProtocol) else "qsp"
    if method == "qsp":
        requires(operator, BlockEncodingProtocol, path="QSP.operator")
        if not callable(qsp):
            raise ValidationError(
                "QSP 路径需要注入实际 qsp(BE,time) 实现；当前库没有通用 QSP-HamSim 内核"
            )
        encoded = as_block_encoding(operator)
        result = qsp(encoded, time)
        require_instance(result, BlockEncoding, "QSP.output")
        if result.width != encoded.width:
            raise ValidationError("QSP 返回的目标宽度不一致")
        return result
    requires(operator, TrotterizableProtocol, path="Trotter.operator")
    positive_integer(steps, "Trotter.steps")
    terms = tuple(operator.trotter_list())
    if not terms:
        raise ValidationError("trotter_list() 需要非空项列表")
    evolutions = []
    for term in terms:
        require_instance(term, TrotterTerm, "Trotter.term")
        # 每项必须提供无后选择的酉演化；一般多项式 BE 不能替代它。
        evolution = term.operator.evolution(term.coefficient * time / steps)
        from pyqecclang.infrastructure.builder import Operation

        require_instance(evolution, Operation, "Trotter.term.evolution")
        if any(r.type.width for r in evolution.module.registers if r.name != "target"):
            raise ValidationError(
                "当前 Trotter 实现要求各项演化仅有 target，其他公开寄存器须为零宽"
            )
        if not any(r.name == "target" for r in evolution.module.registers):
            raise ValidationError("Trotter 项演化需要 target 寄存器")
        evolutions.append(evolution)
    widths = {
        next(r.type.width for r in op.module.registers if r.name == "target") for op in evolutions
    }
    if len(widths) != 1:
        raise ValidationError("Trotter 项演化的目标宽度不一致")
    width = widths.pop()
    b = Builder(
        _name("trotter_protocol", *evolutions, steps),
        {"target": Bits(width), "signal": Bits(0)},
        resources_for(*[(f"term{i}", op) for i, op in enumerate(evolutions)]),
        attributes={
            "algorithm": "trotter_hamiltonian_protocol",
            "steps": steps,
            "correctness": "product formula approximation pending",
        },
    )
    with b.repeat(steps):
        for i, op in enumerate(evolutions):
            arguments = {
                r.name: b["target"] if r.name == "target" else b["signal"]
                for r in op.module.registers
            }
            invoke(b, op, f"term{i}", **arguments)
    return block_encoding(b.finish(), alpha=1.0)
