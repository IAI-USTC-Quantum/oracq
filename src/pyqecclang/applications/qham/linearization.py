"按 HAM 阶数和有序张量字惰性生成 QCL；不物化空间矩阵。"

from __future__ import annotations

import json
import math
from collections.abc import Iterator, Sequence
from dataclasses import dataclass

from pyqecclang.applications.qham.pde import PolynomialPDE
from pyqecclang.infrastructure.ir import ValidationError


def compositions(total: int, length: int) -> Iterator[tuple[int, ...]]:
    """按首分量递增的顺序，惰性枚举和为 ``total`` 的 ``length`` 个非负整数的全部有序拆分。

    Args:
        total: 目标和，为非负整数。
        length: 拆分的元数。

    Returns:
        Iterator[tuple[int, ...]]: 惰性生成器，逐个产出而不物化全部拆分。

    Yields:
        tuple[int, ...]: 长度为 ``length`` 且元素之和等于 ``total`` 的有序元组；
            ``length`` 为 0 时仅当 ``total`` 为 0 产出空元组。
    """
    if length == 0:
        if total == 0:
            yield ()
    elif length == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for rest in compositions(total - first, length - 1):
                yield (first, *rest)


@dataclass(frozen=True)
class HomotopyWeight:
    """QCL 线性边上随同伦参数 ``eta`` 变化的系数。

    Attributes:
        kind: 权重类型，取 ``one``（恒为 1）、``correction`` 或 ``physical``。
        power: ``correction`` 与 ``physical`` 类型所用的幂次。
    """

    kind: str = "one"
    power: int = 0

    def evaluate(self, eta: complex) -> complex:
        """计算权重在同伦参数 ``eta`` 处的取值。

        Args:
            eta: 同伦参数值。

        Returns:
            ``one`` 恒为 1，``correction`` 为 ``-eta*(1+eta)**power``，
            ``physical`` 为 ``1-(1+eta)**power``。

        Raises:
            ValidationError: ``kind`` 不是已支持的权重类型。
        """
        if self.kind == "one":
            return 1.0
        if self.kind == "correction":
            return -eta * (1 + eta) ** self.power
        if self.kind == "physical":
            return 1 - (1 + eta) ** self.power
        raise ValidationError("未知同伦权重")

    def formula(self) -> str:
        """返回权重的公式文本，同伦参数记作 ``eta``。

        Returns:
            str: 与 ``evaluate`` 对应的中缀表达式，如 ``-eta*(1+eta)^3``。
        """
        if self.kind == "one":
            return "1"
        if self.kind == "correction":
            return f"-eta*(1+eta)^{self.power}"
        return f"1-(1+eta)^{self.power}"


@dataclass(frozen=True)
class Block:
    """QCL 有限闭包中的一个块变量。

    Attributes:
        kind: 块类型；``physical`` 为物理输出块（各 HAM 阶分量之和），
            ``tensor`` 为独立坐标上的张量块。
        orders: 张量字，即各因子的 HAM 阶数组成的有序元组；物理块恒为
            空元组，``tensor`` 且为空元组时表示强迫齐次化引入的常量块。
    """

    kind: str
    orders: tuple[int, ...] = ()

    @property
    def rank(self) -> int:
        """块的张量因子个数；物理块按 1 计，空字常量块为 0。"""
        return 1 if self.kind == "physical" else len(self.orders)

    @property
    def label(self) -> str:
        """返回块的可读标签。

        物理块为 ``u_sum``，空字常量块为 ``one``，其余张量块为 ``Y_``
        后接下划线连接的各阶数。
        """
        if self.kind == "physical":
            return "u_sum"
        if not self.orders:
            return "one"
        return "Y_" + "_".join(map(str, self.orders))


@dataclass(frozen=True)
class Coupling:
    """QCL 闭包中的一条线性边，把列块的贡献送入行块的方程。

    Attributes:
        row: 边的目标块（所属的方程行）。
        column: 边的源块（输入变量）。
        operator: 作用的多线性端口名，如 ``L``、``F`` 或 ``B_i``。
        position: 端口作用（替换）行块的第几个张量因子，从 0 计。
        arity: 端口的输入元数。
        weight: 该边的同伦权重，默认恒为 1。
    """

    row: Block
    column: Block
    operator: str
    position: int
    arity: int
    weight: HomotopyWeight = HomotopyWeight()


@dataclass(frozen=True)
class QHAMPlan:
    """多项式 PDE 在 HAM 截断阶下的 QCL 有限闭包计划。

    闭包由权重判据 ``p*sum(orders)+len(orders) <= p*m+1`` 限定，其中
    ``p = max(1, D-1)``、``D`` 为 PDE 多项式次数、``m`` 为截断阶。块集合
    与各行规则均惰性枚举，不物化闭包矩阵。

    Attributes:
        pde: 构造时经过验证的 ``PolynomialPDE``。
        order: HAM 截断阶 ``m``，须为非负整数。
        version: 序列化格式版本，当前为 ``0.1``。

    Raises:
        ValidationError: ``pde`` 未通过验证，或 ``order`` 不是非负整数、
            ``version`` 不是 ``0.1``。
    """

    pde: PolynomialPDE
    order: int
    version: str = "0.1"

    def __post_init__(self) -> None:
        """校验 PDE、截断阶与版本号，非法时抛出 ``ValidationError``。"""
        self.pde.validate()
        if self.version != "0.1" or type(self.order) is not int or self.order < 0:
            raise ValidationError("QHAM 截断阶必须为非负整数")

    @property
    def grade(self) -> int:
        """闭包权重判据的次数缩放 ``p = max(1, D-1)``，``D`` 为 PDE 多项式次数。"""
        return max(1, self.pde.degree - 1)

    @property
    def max_rank(self) -> int:
        """闭包权重上限 ``p*m+1``，即张量块允许的最大张量秩。"""
        return self.grade * self.order + 1

    @property
    def has_forcing(self) -> bool:
        """PDE 是否含常量强迫项（即存在端口 ``F``，闭包需要空字常量块）。"""
        return any(port.name == "F" for port in self.pde.ports)

    def weight(self, word: Sequence[int]) -> int:
        """计算张量字的闭包权重 ``p*sum(word)+len(word)``。

        Args:
            word: 各因子的 HAM 阶数组成的元组。

        Returns:
            int: 该张量字的闭包权重。
        """
        return self.grade * sum(word) + len(word)

    def contains(self, block: Block) -> bool:
        """判断块是否属于该计划的有限闭包。

        物理块须不带阶数；空字张量块仅当存在强迫；其余张量块要求各阶数
        均为非负整数且闭包权重不超过上限；其他块类型一律不在闭包内。

        Args:
            block: 待检查的块。

        Returns:
            bool: 块属于闭包时为 True。
        """
        if block.kind == "physical":
            return block.orders == ()
        if block.kind != "tensor":
            return False
        if not block.orders:
            return self.has_forcing
        return (
            all(type(a) is int and a >= 0 for a in block.orders)
            and self.weight(block.orders) <= self.max_rank
        )

    def rank_count(self, rank: int) -> int:
        """统计给定张量秩的闭包张量块个数。

        Args:
            rank: 张量因子个数，从 1 起计。

        Returns:
            int: 该秩下满足权重约束的张量字数目；秩小于 1 或超过上限时为 0。
        """
        limit = (self.max_rank - rank) // self.grade
        return math.comb(limit + rank, rank) if limit >= 0 and rank >= 1 else 0

    @property
    def block_count(self) -> int:
        """闭包块总数：物理块、可选空字常量块与各秩张量块个数之和。"""
        return (
            1 + int(self.has_forcing) + sum(self.rank_count(k) for k in range(1, self.max_rank + 1))
        )

    def blocks(self, *, max_blocks: int | None = None) -> Iterator[Block]:
        """按规范顺序惰性枚举闭包中的全部块。

        物理块在前，张量块按秩升序、同秩内按字和升序（同字和按 ``compositions``
        的枚举序）排列，空字常量块（若有强迫）在最后。

        Args:
            max_blocks: 显式枚举的块数预算；块总数超过该值时立即报错，
                提示改用惰性查询。

        Returns:
            Iterator[Block]: 惰性生成器，逐个产出闭包块而不物化完整块列表。

        Yields:
            Block: 闭包中的块。

        Raises:
            ValidationError: 指定了 ``max_blocks`` 且闭包块数超过预算。
        """
        if max_blocks is not None and self.block_count > max_blocks:
            raise ValidationError("显式块生成超过预算；QHAMPlan 仍可惰性查询")
        yield Block("physical")
        for rank in range(1, self.max_rank + 1):
            limit = (self.max_rank - rank) // self.grade
            for total in range(limit + 1):
                for word in compositions(total, rank):
                    yield Block("tensor", word)
        if self.has_forcing:
            yield Block("tensor")

    def row_terms(self, row: Block) -> tuple[Coupling, ...]:
        """枚举行块在闭包内的全部线性边。

        覆盖三类边：同字的线性边 ``L``；强迫插入边 ``F``（物理行来自空字
        常量块，张量行来自删去该零阶因子后的字）；把高阶因子替换为低阶
        多元组并施加非线性端口的边（携带 ``physical`` 或 ``correction``
        同伦权重）。空字常量块的行没有边。

        Args:
            row: 行块，须属于闭包。

        Returns:
            tuple[Coupling, ...]: 该行的全部边，每条边的 ``row`` 字段即入参。

        Raises:
            ValidationError: 行块不属于有限闭包。
            AssertionError: 生成的列块越出闭包，即闭包构造不变量被破坏。
        """
        if not self.contains(row):
            raise ValidationError("QCL 行不属于有限闭包")
        if row.kind == "tensor" and not row.orders:
            return ()
        ports = self.pde.ports
        linear = any(p.name == "L" for p in ports)
        nonlinear = [p for p in ports if p.arity >= 2]
        result: list[Coupling] = []
        if row.kind == "physical":
            if linear:
                result.append(Coupling(row, row, "L", 0, 1))
            if self.has_forcing:
                result.append(Coupling(row, Block("tensor"), "F", 0, 0))
            for total in range(self.order):
                for port in nonlinear:
                    for word in compositions(total, port.arity):
                        result.append(
                            Coupling(
                                row,
                                Block("tensor", word),
                                port.name,
                                0,
                                port.arity,
                                HomotopyWeight("physical", self.order - total),
                            )
                        )
        else:
            for position, order in enumerate(row.orders):
                if linear:
                    result.append(Coupling(row, row, "L", position, 1))
                if order == 0:
                    if self.has_forcing:
                        word = row.orders[:position] + row.orders[position + 1 :]
                        result.append(Coupling(row, Block("tensor", word), "F", position, 0))
                else:
                    for total in range(order):
                        for port in nonlinear:
                            for word in compositions(total, port.arity):
                                source = row.orders[:position] + word + row.orders[position + 1 :]
                                result.append(
                                    Coupling(
                                        row,
                                        Block("tensor", source),
                                        port.name,
                                        position,
                                        port.arity,
                                        HomotopyWeight("correction", order - 1 - total),
                                    )
                                )
        if any(not self.contains(term.column) for term in result):
            raise AssertionError("QCL 闭包构造错误")
        return tuple(result)

    def raw_dimension(self, dimension: int) -> int:
        """计算闭包线性系统的总维数。

        Args:
            dimension: 单块局部状态空间的维数。

        Returns:
            int: 物理块与各秩张量块的局部维数（``dimension**rank``）之和，
            再加空字常量块的 1 个坐标（若有强迫）。
        """
        return (
            dimension
            + sum(self.rank_count(k) * dimension**k for k in range(1, self.max_rank + 1))
            + int(self.has_forcing)
        )

    def offset(self, block: Block, dimension: int) -> int:
        """返回块内局部坐标 0 在闭包系统中的全局起始索引。

        块按 ``blocks`` 的规范顺序展开：物理块在最前，各秩张量块占
        ``dimension**rank`` 个连续坐标，空字常量块占最后 1 个坐标。

        Args:
            block: 闭包内的块。
            dimension: 单块局部状态空间的维数。

        Returns:
            int: 该块的起始全局索引。

        Raises:
            ValidationError: 块不属于闭包。
        """
        if not self.contains(block):
            raise ValidationError("未知 QCL 块")
        if block.kind == "physical":
            return 0
        if not block.orders:
            return self.raw_dimension(dimension) - 1
        rank, total = len(block.orders), sum(block.orders)
        offset = dimension + sum(self.rank_count(k) * dimension**k for k in range(1, rank))
        before = math.comb(total + rank - 1, rank) if total else 0
        remainder = total
        for pos, value in enumerate(block.orders[:-1]):
            slots = rank - pos - 1
            for first in range(value):
                before += math.comb(remainder - first + slots - 1, slots - 1)
            remainder -= value
        return offset + before * dimension**rank

    def locate(self, index: int, dimension: int) -> tuple[Block, int]:
        """把全局索引分解为所属块与块内局部索引，为 ``offset`` 的逆。

        Args:
            index: 闭包系统内的全局索引。
            dimension: 单块局部状态空间的维数。

        Returns:
            tuple[Block, int]: 索引所在的块及块内局部索引。

        Raises:
            ValidationError: 索引越出闭包系统总维数。
        """
        if not 0 <= index < self.raw_dimension(dimension):
            raise ValidationError("QCL 原始索引越界")
        if index < dimension:
            return Block("physical"), index
        index -= dimension
        for rank in range(1, self.max_rank + 1):
            size = dimension**rank
            group = self.rank_count(rank) * size
            if index >= group:
                index -= group
                continue
            ordinal, local = divmod(index, size)
            total = 0
            while ordinal >= math.comb(total + rank - 1, rank - 1):
                ordinal -= math.comb(total + rank - 1, rank - 1)
                total += 1
            word: list[int] = []
            remaining = total
            for pos in range(rank - 1):
                slots = rank - pos - 1
                first = 0
                while ordinal >= math.comb(remaining - first + slots - 1, slots - 1):
                    ordinal -= math.comb(remaining - first + slots - 1, slots - 1)
                    first += 1
                word.append(first)
                remaining -= first
            word.append(remaining)
            return Block("tensor", tuple(word)), local
        return Block("tensor"), 0

    def summary(self) -> dict[str, object]:
        """汇总计划的关键量，供报告与诊断使用。

        Returns:
            dict: 含版本、截断阶、PDE 次数、最大张量秩、块总数、是否含
            强迫、闭包判据描述与各端口元数的字典。
        """
        return {
            "version": self.version,
            "order": self.order,
            "pde_degree": self.pde.degree,
            "max_tensor_rank": self.max_rank,
            "block_count": self.block_count,
            "homogeneous_constant": self.has_forcing,
            "closure": "(D-1)*sum(orders)+rank <= (D-1)*m+1; D<=1 uses grade 1",
            "operators": [
                {"name": p.name, "input_rank": p.arity, "output_rank": 1} for p in self.pde.ports
            ],
        }

    def dumps(self) -> str:
        """把计划序列化为规范 JSON 文本。

        Returns:
            str: 含 ``version``、``order`` 与 ``pde`` 字段的 JSON 文本，
            键排序、缩进为 2 且以换行结尾。
        """
        return (
            json.dumps(
                {"version": self.version, "order": self.order, "pde": json.loads(self.pde.dumps())},
                ensure_ascii=False,
                sort_keys=True,
                indent=2,
            )
            + "\n"
        )

    @classmethod
    def loads(cls, text: str) -> QHAMPlan:
        """从 ``dumps`` 输出的 JSON 文本重建计划。

        Args:
            text: JSON 文本。

        Returns:
            QHAMPlan: 重建后的计划，构造时会再次执行验证。

        Raises:
            ValidationError: JSON 结构、字段集合或数据无效。
        """
        try:
            raw = json.loads(text)
            if set(raw) != {"version", "order", "pde"}:
                raise ValidationError("未知 QCL plan 字段")
            return cls(PolynomialPDE.loads(json.dumps(raw["pde"])), raw["order"], raw["version"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValidationError("无效 QCL plan：" + str(exc)) from exc
