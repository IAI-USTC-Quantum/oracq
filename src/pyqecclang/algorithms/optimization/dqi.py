"""DQI（Decoded Quantum Interferometry，解码量子干涉优化）的 GF(2) max-XORSAT 生成器。

实现 Jordan et al. 2024（arXiv:2408.08292，图 4）的线路骨架：在 m 比特
error 寄存器制备权重 l 的 Dicke 态，施加右端项相位 ``(-1)**(v.y)``，把 ``B^T y``
可逆计算进 n 比特 syndrome 寄存器，再用可逆经典译码器把 error 寄存器卸载回
``|0>``，最后对 syndrome 寄存器做 Hadamard 变换并测量。译码器是输入模型的一
部分，以开放声明（abstract）加分批绑定（bind）的三层范式接入。

当前只支持 GF(2)：GF(q) 情形需要 q 元离散 Fourier 变换与广义 Dicke 态，
寄存器还需按 ``log2(q)`` 分子寄存器组织，留作扩展。测量后保留 error 寄存器
为零的分支（后选）由调用方完成。
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from math import comb

from pyqecclang.algorithms.input_model.contracts import (
    OracleView,
    fail,
    positive_integer,
    require_instance,
    validate_signature,
)
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import (
    StatePreparation,
    annotate,
    declare,
    gate_state_prep,
    invoke,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder, Operation
from pyqecclang.infrastructure.ir import Bits, ValidationError

# 显式幅度/穷举译码只服务小实例见证，超出该规模应换成专用制备线路或高效译码器。
_EXPLICIT_LIMIT = 16


@dataclass(frozen=True)
class XorSatInstance:
    """GF(2) max-XORSAT 实例：约束组 ``Bx = v`` 的稀疏行表示。

    rows 每行一个约束，内容为参与该约束的变量下标，行内不重复；rhs 为右端项 ``v``，
    长度与行数相同，取值 0 或 1；num_variables 为变量数 n，也是 syndrome 寄存器的位宽。
    """

    rows: tuple[tuple[int, ...], ...]
    rhs: tuple[int, ...]
    num_variables: int

    def __post_init__(self) -> None:
        """规范化并校验约束行、右端项与变量下标的构造期约束。"""
        positive_integer(self.num_variables, "XorSatInstance.num_variables", maximum=64)
        rows = tuple(tuple(row) for row in self.rows)
        rhs = tuple(self.rhs)
        if not rows or len(rows) != len(rhs):
            fail(
                "CONFIG_VALUE",
                "XorSatInstance.rows",
                "非空且与 rhs 等长",
                (len(rows), len(rhs)),
                "约束行必须非空且与右端项等长",
            )
        for row in rows:
            if not row or len(set(row)) != len(row):
                fail(
                    "CONFIG_VALUE",
                    "XorSatInstance.rows",
                    "非空且无重复下标",
                    row,
                    "约束行必须非空且行内变量下标不重复",
                )
            for index in row:
                positive_integer(
                    index,
                    "XorSatInstance.index",
                    minimum=0,
                    maximum=self.num_variables - 1,
                )
        if any(type(bit) is not int or bit not in (0, 1) for bit in rhs):
            fail(
                "CONFIG_VALUE",
                "XorSatInstance.rhs",
                "0 或 1",
                rhs,
                "右端项必须取 0 或 1，不能是 bool",
            )
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "rhs", rhs)

    @property
    def num_constraints(self) -> int:
        """约束数 m，即 error 寄存器的位宽。"""
        return len(self.rows)

    def satisfied_count(self, assignment: int) -> int:
        """统计给定整数赋值满足的约束数。

        Args:
            assignment: 赋值的整数编码，取 0..2^变量数−1，第 j 位是变量 j 的取值。

        Returns:
            int: 该赋值满足的约束条数。
        """
        positive_integer(
            assignment,
            "XorSatInstance.assignment",
            minimum=0,
            maximum=(1 << self.num_variables) - 1,
        )
        return sum(
            (sum((assignment >> j) & 1 for j in row) & 1) == bit
            for row, bit in zip(self.rows, self.rhs, strict=True)
        )


@dataclass(frozen=True)
class DecoderOracle(OracleView):
    """综合征译码器的可逆形式，是 DQI 的关键开放输入。

    语义为 ``|syndrome, error> -> |syndrome, error XOR D(syndrome)>``，其中
    ``D`` 是经典译码函数（例如 belief propagation 的可逆实现）。无法译码的
    综合征可以映射到任意错误模式，对应分支在后选中被淘汰。
    """

    oracle_kind = "reversible_function"
    operation: Operation

    def decoder(self) -> DecoderOracle:
        """译码器角色访问器，返回自身；与其他 ``OracleView`` 的角色方法一致。

        Returns:
            DecoderOracle: 自身引用，保持角色访问器接口一致。
        """
        return self

    def __post_init__(self) -> None:
        """校验包装操作具有 syndrome 与 error 签名。"""
        validate_signature(self.operation, ("syndrome", "error"), "DecoderOracle")

    @property
    def syndrome_width(self) -> int:
        """syndrome 寄存器位宽，直接读取 RIR 寄存器签名。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "syndrome")

    @property
    def error_width(self) -> int:
        """error 寄存器位宽，直接读取 RIR 寄存器签名。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "error")


def abstract_decoder(name: str, syndrome_width: int, error_width: int) -> DecoderOracle:
    """声明译码器槽位；高效经典译码算法经 bind 分批绑定。

    Args:
        name: 译码器槽位的声明模块名。
        syndrome_width: syndrome 寄存器位宽，取 1..64。
        error_width: error 寄存器位宽，取 1..64。

    Returns:
        DecoderOracle: 体为空、由 bind 延迟绑定实现的译码器槽位句柄。
    """
    positive_integer(syndrome_width, "abstract_decoder.syndrome_width", maximum=64)
    positive_integer(error_width, "abstract_decoder.error_width", maximum=64)
    return DecoderOracle(
        declare(
            name,
            {"syndrome": Bits(syndrome_width), "error": Bits(error_width)},
            paradigm="reversible_function",
            attributes={"decoder_role": "syndrome_decode"},
        )
    )


def table_decoder(
    syndrome_width: int, error_width: int, table: Mapping[int, int], *, name: str | None = None
) -> DecoderOracle:
    """用显式查询表实现译码器的 gate 见证；只译码表中列出的综合征。

    Args:
        syndrome_width: syndrome 寄存器位宽，即实例的变量数 n。
        error_width: error 寄存器位宽，即实例的约束数 m。
        table: 综合征整数值到错误模式整数值的映射。
        name: 可选模块名。

    Returns:
        DecoderOracle: 未列出的综合征映射到零错误模式。"""
    positive_integer(syndrome_width, "table_decoder.syndrome_width", maximum=64)
    positive_integer(error_width, "table_decoder.error_width", maximum=64)
    items = tuple(sorted(table.items()))
    b = Builder(
        name or _name("decoder_table", syndrome_width, error_width, items),
        {"syndrome": Bits(syndrome_width), "error": Bits(error_width)},
    )
    for syndrome, error in items:
        if (
            type(syndrome) is not int
            or type(error) is not int
            or not (0 <= syndrome < 1 << syndrome_width and 0 <= error < 1 << error_width)
        ):
            raise ValidationError("译码查询表的综合征或错误模式越界")
        if error:
            with b.control(b["syndrome"], syndrome):
                for bit in range(error_width):
                    if (error >> bit) & 1:
                        b.x(b["error"][bit])
    return DecoderOracle(
        annotate(b.finish(), "reversible_function", implementation="gate_truth_table")
    )


def bruteforce_decoder(
    instance: XorSatInstance, *, max_weight: int | None = None, name: str | None = None
) -> DecoderOracle:
    """穷举最小权重译码器：为每个综合征给出权重不超过 max_weight 的最轻错误。

    枚举全部 2**m 个错误模式，只接受 m <= 16 的小实例；大实例应绑定高效
    译码器的可逆实现。等权重并列时保留枚举序最小的错误模式。

    Args:
        instance: 被译码的 XOR 实例，约束数不超过 16。
        max_weight: 错误模式的权重上限，取 0..约束数；缺省为约束数。
        name: 可选的生成模块名。

    Returns:
        DecoderOracle: 每个综合征映射到最轻错误模式的显式查表译码器。
    """
    require_instance(instance, XorSatInstance, "bruteforce_decoder.instance")
    m = instance.num_constraints
    if m > _EXPLICIT_LIMIT:
        fail(
            "CONFIG_VALUE",
            "bruteforce_decoder.m",
            f"<= {_EXPLICIT_LIMIT}",
            m,
            "穷举译码器只接受小规模实例",
        )
    if max_weight is None:
        max_weight = m
    positive_integer(max_weight, "bruteforce_decoder.max_weight", minimum=0, maximum=m)
    row_masks = [sum(1 << j for j in row) for row in instance.rows]
    best: dict[int, int] = {}
    for error in range(1 << m):
        weight = error.bit_count()
        if weight > max_weight:
            continue
        syndrome = 0
        for i, mask in enumerate(row_masks):
            if (error >> i) & 1:
                syndrome ^= mask
        if syndrome not in best or weight < best[syndrome].bit_count():
            best[syndrome] = error
    return table_decoder(
        instance.num_variables,
        m,
        best,
        name=name
        or _name(
            "decoder_bruteforce",
            instance.rows,
            instance.rhs,
            instance.num_variables,
            max_weight,
        ),
    )


def dicke_state(m: int, weight: int) -> StatePreparation:
    """制备 m 比特、权重 weight 的 Dicke 态 ``|D_l^m>``。

    采用显式幅度的多路旋转构造，只接受 m <= 16 的小规模；大规模制备需要
    专门的 Dicke 态线路（如 Bartschi-Eidenbenz 的 O(l*m) 构造），可作为
    state_prep_isometry 开放声明另行接入。

    Args:
        m: 目标位数，范围 1..16。
        weight: Hamming 权重，范围 0..m。

    Returns:
        StatePreparation: 零输入制备，target 位宽为 m。"""
    positive_integer(m, "dicke_state.m", maximum=_EXPLICIT_LIMIT)
    positive_integer(weight, "dicke_state.weight", minimum=0, maximum=m)
    amplitude = 1 / math.sqrt(comb(m, weight))
    amplitudes = [amplitude if value.bit_count() == weight else 0 for value in range(1 << m)]
    return gate_state_prep(amplitudes, name=f"dicke_{m}_{weight}")


def dqi(instance: XorSatInstance, decoder: DecoderOracle, *, weight: int) -> Operation:
    """组装 DQI 主线路（GF(2)，论文图 4 的单权重版本）。

    Args:
        instance: XorSatInstance，约束数 m 与变量数 n 决定寄存器位宽。
        decoder: DecoderOracle，接口宽度必须与实例的 n/m 一致；可以是抽象声明。
        weight: Dicke 态权重 l，范围 0..m；译码半径需覆盖该权重才有意义。

    Returns:
        Operation: 寄存器为 error(m) 与 syndrome(n)。从全零初态出发，在译码
        成功的分支上 error 回到 ``|0>``，syndrome 经 Hadamard 后以正比于
        ``K_l(u(x))**2`` 的概率测得赋值 x，其中 u(x) 是未满足约束数、``K_l``
        是 Krawtchouk 多项式。后选 error 为零由调用方完成。"""
    require_instance(instance, XorSatInstance, "dqi.instance")
    require_instance(decoder, DecoderOracle, "dqi.decoder")
    m, n = instance.num_constraints, instance.num_variables
    positive_integer(weight, "dqi.weight", minimum=0, maximum=m)
    if decoder.syndrome_width != n or decoder.error_width != m:
        fail(
            "INPUT_WIDTH",
            "dqi.decoder",
            (n, m),
            (decoder.syndrome_width, decoder.error_width),
            "译码器位宽必须与实例的 n/m 一致",
        )
    preparation = dicke_state(m, weight)
    b = Builder(
        _name("dqi", instance, decoder.operation, weight),
        {"error": Bits(m), "syndrome": Bits(n)},
        resources_for(("prep", preparation.operation), ("decoder", decoder.operation)),
        attributes={
            "algorithm": "dqi",
            "field": "GF(2)",
            "num_constraints": m,
            "num_variables": n,
            "dicke_weight": weight,
            "readout_register": "syndrome",
            "postselection": "error_zero",
        },
    )
    invoke(b, preparation.operation, "prep", target=b["error"], work=b["syndrome"][:0])
    for i, bit in enumerate(instance.rhs):
        if bit:
            b.z(b["error"][i])
    for i, row in enumerate(instance.rows):
        for j in row:
            b.xor(b["error"][i], b["syndrome"][j])
    invoke(b, decoder.operation, "decoder", syndrome=b["syndrome"], error=b["error"])
    b.h(b["syndrome"])
    return b.finish()
