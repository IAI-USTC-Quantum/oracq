"LCU 的 PREPARE–SELECT 标准分解（Low & Chuang 2019；Babbush et al. 2018 alias 采样）。"

from __future__ import annotations

import cmath
import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from pyqecclang.algorithms.common.arithmetic import FixedFormat, fixed_arithmetic
from pyqecclang.algorithms.common.hamiltonian import PauliHamiltonian
from pyqecclang.algorithms.input_model.block_encoding import pauli_word
from pyqecclang.algorithms.input_model.contracts import positive_integer, require_instance
from pyqecclang.algorithms.input_model.interfaces import (
    StatePreparationProtocol,
    as_state_preparation,
)
from pyqecclang.algorithms.input_model.operators import BlockEncoding, _name, scale
from pyqecclang.algorithms.input_model.oracles import (
    StatePreparation,
    XorDatabase,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    qram_database,
    qram_state_angles,
    qram_state_prep,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, ValidationError


def _normalized(
    coefficients: Iterable[complex],
) -> tuple[tuple[complex, ...], float, int, list[float]]:
    """校验系数并给出 (系数, alpha=l1 范数, selector 位宽, ``√|c|/√α`` 振幅)。"""
    values = tuple(complex(c) for c in coefficients)
    if len(values) < 2:
        raise ValidationError("PREPARE 至少需要两个系数；单项请直接用 scale")
    if not all(math.isfinite(c.real) and math.isfinite(c.imag) for c in values):
        raise ValidationError("LCU 系数必须有限")
    alpha = sum(abs(c) for c in values)
    if not alpha:
        raise ValidationError("LCU 系数不能全为零")
    width = (len(values) - 1).bit_length()
    amplitudes = [math.sqrt(abs(c) / alpha) for c in values]
    amplitudes += [0.0] * ((1 << width) - len(amplitudes))
    return values, alpha, width, amplitudes


def _pauli_terms(
    terms: Iterable[tuple[complex, str]] | PauliHamiltonian,
) -> tuple[tuple[complex, str], ...]:
    """接受 (系数, Pauli 字) 序列或 PauliHamiltonian；保留零系数以对齐索引。"""
    from pyqecclang.algorithms.common.hamiltonian import PauliHamiltonian

    if isinstance(terms, PauliHamiltonian):
        terms = terms.terms
    terms = tuple((complex(c), str(w)) for c, w in terms)
    if not terms:
        raise ValidationError("LCU 至少需要一个 Pauli 项")
    for coefficient, word in terms:
        if not (math.isfinite(coefficient.real) and math.isfinite(coefficient.imag)):
            raise ValidationError("Pauli 系数必须有限")
        if not word or any(letter not in "IXYZ" for letter in word):
            raise ValidationError("Pauli 字只允许 I/X/Y/Z")
    if len({len(word) for _, word in terms}) != 1:
        raise ValidationError("Pauli 项宽度不一致")
    return terms


def _prepare_attributes(
    values: Sequence[complex], alpha: float, width: int
) -> dict[str, int | float]:
    """汇总 PREPARE 写入模块属性的项数、alpha 与 selector 位宽。"""
    return {
        "prepare_terms": len(values),
        "prepare_alpha": alpha,
        "selector_width": width,
    }


def abstract_prepare(
    coefficients: Iterable[complex], *, work_width: int = 0, name: str | None = None
) -> StatePreparation:
    """PREPARE 的开放声明：系数写入声明属性，体为空，由 bind 延迟绑定实现。

    work_width 必须与后续绑定实现的 work 宽度一致（gate 为 0，QRAM 为
    selector 位宽加角度位宽），与 catalog 中 abstract_state_prep 的用法一致。

    Args:
        coefficients: LCU 复系数序列，至少两项且取值有限，不能全为零。
        work_width: work 寄存器位宽，取 0..64；须与后续绑定实现的 work 宽度一致。
        name: 声明模块名；缺省按系数与位宽自动生成。

    Returns:
        StatePreparation: 包装开放声明 PREPARE 模块的状态制备句柄。
    """
    values, alpha, width, _ = _normalized(coefficients)
    positive_integer(work_width, "abstract_prepare.work_width", minimum=0)
    return StatePreparation(
        declare(
            name or _name("prepare_abstract", values, work_width),
            {"target": Bits(width), "work": Bits(work_width)},
            paradigm="state_prep_isometry",
            attributes={
                "zero_input": True,
                "clean_work": True,
                **_prepare_attributes(values, alpha, width),
            },
        )
    )


def gate_prepare(
    coefficients: Iterable[complex], *, name: str | None = None
) -> StatePreparation:
    """门级 PREPARE：振幅 ∝ ``√|c_i|`` 的多路复用 Ry；系数相位按仓库惯例进 SELECT。

    Args:
        coefficients: LCU 复系数序列，至少两项且取值有限，不能全为零。
        name: 生成的状态制备模块名；缺省按系数自动生成。

    Returns:
        StatePreparation: 门级实现的状态制备操作句柄。
    """
    values, alpha, width, amplitudes = _normalized(coefficients)
    prep = gate_state_prep(amplitudes, name=name)
    return StatePreparation(
        annotate(
            prep.operation,
            "state_prep_isometry",
            zero_input=True,
            clean_work=True,
            **_prepare_attributes(values, alpha, width),
        )
    )


@dataclass(frozen=True)
class QramPreparation:
    """QRAM 版 PREPARE 句柄；角度表在 memory 中提供，不进入 IR。"""

    preparation: StatePreparation
    memory: dict[str, dict[int, int]]

    def state_preparation(self) -> StatePreparation:
        """返回句柄内包装的 ``StatePreparation`` 操作。

        Returns:
            StatePreparation: 句柄内包装的状态制备操作。
        """
        return self.preparation


def qram_prepare(coefficients: Iterable[complex], *, angle_width: int = 8) -> QramPreparation:
    """QRAM 资源版 PREPARE：旋转角度表由调用方作为 QRAM 数据绑定。

    Args:
        coefficients: LCU 复系数序列，至少两项且取值有限，不能全为零。
        angle_width: 旋转角的量化位宽，取 2..32。

    Returns:
        QramPreparation: 含 PREPARE 操作与默认角度表 QRAM 数据的句柄。
    """
    values, alpha, width, amplitudes = _normalized(coefficients)
    positive_integer(angle_width, "qram_prepare.angle_width", minimum=2, maximum=32)
    prep = qram_state_prep(width, angle_width)
    memory = {"angles": qram_state_angles(amplitudes, angle_width)}
    operation = annotate(
        prep.operation,
        "state_prep_isometry",
        zero_input=True,
        clean_work=True,
        **_prepare_attributes(values, alpha, width),
    )
    return QramPreparation(StatePreparation(operation), memory)


@dataclass(frozen=True)
class AliasTable:
    """Babbush et al. alias 采样的经典预处理结果；字打包为 ``keep | (alt << precision)``。"""

    probabilities: tuple[float, ...]
    keep: tuple[float, ...]
    alt: tuple[int, ...]
    quantized: tuple[int, ...]
    precision: int
    table: dict[int, int]

    def distribution(self) -> tuple[float, ...]:
        """量化 keep 与均匀抽取下的经典采样分布。

        Returns:
            tuple[float, ...]: 各槽位的采样概率，长度等于槽总数，总和为一。
        """
        scale = 1 << self.precision
        size = len(self.keep)
        result = [0] * size
        for i in range(size):
            result[i] += self.quantized[i]
            result[self.alt[i]] += scale - self.quantized[i]
        return tuple(v / (size * scale) for v in result)


def alias_table(coefficients: Iterable[complex], *, precision: int = 8) -> AliasTable:
    """Vose alias 预处理：padded 到 2^selector 后按均值分裂 keep/alt；零概率槽也可工作。

    Args:
        coefficients: LCU 复系数序列，至少两项且取值有限，不能全为零。
        precision: keep 概率的量化位宽，取 2..32；selector 位宽与之的和不得超过 64。

    Returns:
        AliasTable: 含概率、keep/alt 表、量化字与打包数据表的预处理结果。
    """
    values, alpha, width, _ = _normalized(coefficients)
    positive_integer(precision, "alias_table.precision", minimum=2, maximum=32)
    size = 1 << width
    if width + precision > 64:
        raise ValidationError("alias 数据字超出 QRAM 位宽上限")
    probabilities = [abs(c) / alpha for c in values] + [0.0] * (size - len(values))
    scaled = [p * size for p in probabilities]
    keep, alt = [0.0] * size, list(range(size))
    small = [i for i, v in enumerate(scaled) if v < 1]
    large = [i for i, v in enumerate(scaled) if v >= 1]
    while small and large:
        lo, hi = small.pop(), large.pop()
        keep[lo] = scaled[lo]
        alt[lo] = hi
        scaled[hi] += scaled[lo] - 1
        (small if scaled[hi] < 1 else large).append(hi)
    for i in small + large:
        keep[i] = 1.0
    quantized = [min((1 << precision) - 1, math.floor(k * (1 << precision))) for k in keep]
    table = {i: quantized[i] | (alt[i] << precision) for i in range(size)}
    return AliasTable(
        tuple(probabilities), tuple(keep), tuple(alt), tuple(quantized), precision, table
    )


@dataclass(frozen=True)
class AliasPreparation:
    """alias 采样 PREPARE 句柄；table 是经典表，memory 是默认 QRAM 绑定的数据。"""

    preparation: StatePreparation
    table: AliasTable
    memory: dict[str, dict[int, int]]

    def state_preparation(self) -> StatePreparation:
        """返回句柄内包装的 ``StatePreparation`` 操作。

        Returns:
            StatePreparation: 句柄内包装的状态制备操作。
        """
        return self.preparation


def alias_prepare(
    coefficients: Iterable[complex],
    *,
    precision: int = 8,
    database: XorDatabase | None = None,
) -> AliasPreparation:
    """Babbush et al. alias 采样 PREPARE：均匀态 + (keep,alt) 加载 + 比较器 + 受控交换。

    数据表通过 XorDatabase 接口注入（默认 QRAM 资源，也可传 gate_database），
    不嵌入 IR。work 中的数据/比较位与 selector 纠缠，构成 ``Σ√p_i|i>|junk_i>`` 的
    纯化 junk；块编码 (0,0) 块不受 junk 内积影响，junk 由伴随 PREPARE 复净。
    keep 按 precision 位向下取整量化，与精确分布的总变差不超过 2^selector·2^-precision。

    Args:
        coefficients: LCU 复系数序列，至少两项且取值有限，不能全为零。
        precision: keep 概率的量化位宽，取 2..32。
        database: 承载 (keep|alt) 数据表的 XorDatabase 句柄；缺省为 QRAM 资源，
            地址位宽须等于 selector 位宽、数据位宽等于 precision+selector 位宽。

    Returns:
        AliasPreparation: 含 PREPARE 操作、经典 alias 表与默认 QRAM 数据的句柄。
    """
    values, alpha, width, _ = _normalized(coefficients)
    positive_integer(precision, "alias_prepare.precision", minimum=2, maximum=32)
    table = alias_table(coefficients, precision=precision)
    data_width = precision + width
    if database is None:
        database = qram_database(width, data_width)
        memory = {
            f"db__{r.name}": dict(table.table) for r in database.operation.module.resources
        }
    else:
        require_instance(database, XorDatabase, "alias_prepare.database")
        memory = {}
    if database.address_width != width or database.data_width != data_width:
        raise ValidationError("alias 数据表需要 address=selector、``data=(keep|alt)`` 的宽度")
    compare = fixed_arithmetic("lt", FixedFormat(precision, 0, signed=False))
    b = Builder(
        _name("alias_prepare", values, precision, database.operation),
        {"target": Bits(width), "work": Bits(data_width + precision)},
        resources_for(("db", database.operation), ("lt", compare)),
    )
    data, coin = b["work"][:data_width], b["work"][data_width:]
    keep, alt = data[:precision], data[precision:]
    flag, status = b.local("flag", Bits(1)), b.local("status", Bits(2))
    b.h(b["target"])
    b.h(coin)
    invoke(b, database.operation, "db", address=b["target"], data=data)
    invoke(b, compare, "lt", a=coin, b=keep, out=flag, status=status)
    with b.control(flag, 0):
        b.swap(b["target"], alt)
    # 比较输入未因交换改变，再次调用即可复净 flag；data/coin 作为纯化 junk 留在 work。
    invoke(b, compare, "lt", a=coin, b=keep, out=flag, status=status)
    operation = annotate(
        b.finish(),
        "state_prep_isometry",
        zero_input=True,
        clean_work=False,
        implementation="alias_sampling",
        alias_precision=precision,
        **_prepare_attributes(values, alpha, width),
    )
    return AliasPreparation(StatePreparation(operation), table, memory)


def select_pauli(terms: Iterable[tuple[complex, str]] | PauliHamiltonian) -> Operation:
    """SELECT：selector==i 时对 target 施加第 i 个 Pauli 字，并施加系数相位。

    控制条件按 selector 的二进制值用 RIR Control 原语表达；需要一元迭代
    （unary iteration）时由后端把多比特控制降级实现，生成阶段不展开。

    Args:
        terms: (系数, Pauli 字) 序列或 PauliHamiltonian；至少两项，Pauli 字
            只含 I/X/Y/Z 且等宽。

    Returns:
        Operation: 按 selector 取值施加对应 Pauli 字与系数相位的 SELECT 操作。
    """
    terms = _pauli_terms(terms)
    if len(terms) < 2:
        raise ValidationError("SELECT 至少需要两个 Pauli 项")
    width = len(terms[0][1])
    selector_width = (len(terms) - 1).bit_length()
    b = Builder(
        _name("select_pauli", terms),
        {"selector": Bits(selector_width), "target": Bits(width)},
    )
    for index, (coefficient, word) in enumerate(terms):
        with b.control(b["selector"], index):
            if cmath.phase(coefficient):
                b.global_phase(cmath.phase(coefficient))
            for bit, letter in enumerate(word):
                if letter != "I":
                    b.gate(letter.lower(), b["target"][bit])
    return annotate(
        b.finish(),
        "unitary",
        algorithm="select_pauli",
        select_terms=len(terms),
        selector_width=selector_width,
    )


def lcu_prepare_select(
    terms: Iterable[tuple[complex, str]] | PauliHamiltonian,
    *,
    prepare: StatePreparationProtocol | None = None,
) -> BlockEncoding:
    """标准 PREPARE–SELECT 块编码：(PREPARE†⊗I)·SELECT·(PREPARE⊗I)。

    signal 低 selector_width 位是 selector，高位是 PREPARE 的 work（alias 版为
    纯化 junk）；prepare 缺省为 gate_prepare，也可传 abstract/qram/alias 句柄。
    单项退化为 scale。可直接交给 transforms.qubitization_walk。

    Args:
        terms: (系数, Pauli 字) 序列或 PauliHamiltonian；Pauli 字只含 I/X/Y/Z 且等宽。
        prepare: PREPARE 句柄，缺省为 gate_prepare；须满足 zero_input 契约且
            目标宽度等于 selector 位宽。

    Returns:
        BlockEncoding: 尺度为系数 l1 范数的 (PREPARE†⊗I)·SELECT·(PREPARE⊗I) 块编码。
    """
    terms = _pauli_terms(terms)
    if len(terms) == 1:
        return scale(terms[0][0], pauli_word(terms[0][1]))
    coefficients = tuple(c for c, _ in terms)
    _, alpha, selector_width, _ = _normalized(coefficients)
    prep = gate_prepare(coefficients) if prepare is None else as_state_preparation(prepare)
    if prep.width != selector_width:
        raise ValidationError("PREPARE 目标宽度与 selector 位宽不匹配")
    if dict(prep.operation.module.attributes).get("zero_input") is not True:
        raise ValidationError("PREPARE 需要 zero_input 契约")
    select = select_pauli(terms)
    width = len(terms[0][1])
    junk = prep.work_width
    b = Builder(
        _name("lcu_prepare_select", select, prep.operation),
        {"target": Bits(width), "signal": Bits(selector_width + junk)},
        resources_for(("prep", prep.operation), ("select", select)),
    )
    selector, work = b["signal"][:selector_width], b["signal"][selector_width:]
    invoke(b, prep.operation, "prep", target=selector, work=work)
    invoke(b, select, "select", selector=selector, target=b["target"])
    with b.adjoint():
        invoke(b, prep.operation, "prep", target=selector, work=work)
    return BlockEncoding(
        annotate(
            b.finish(),
            "block_encoding",
            be_alpha=alpha,
            be_form="prepare_select",
            lcu_terms=len(terms),
            selector_width=selector_width,
            prepare_work=junk,
        )
    )
