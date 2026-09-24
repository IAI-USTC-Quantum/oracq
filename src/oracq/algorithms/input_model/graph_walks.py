"""图 oracle input model 与 Szegedy/MNRS 量子行走搜索框架。

图以邻居表 oracle 输入：``|v>|j>|c> -> |v>|j>|c XOR N(v,j)>``，其中出边表统一补齐到
D=2**degree_bits 列（通常用自环补齐）。Szegedy 行走步定义在 (current, peer, index)
二部行走空间上，由邻居叠加态制备的两次反射构成 W = R_B R_A；MNRS 搜索骨架为
setup 初态制备后 Repeat{ walk; marked 相位翻转 }。

参考文献：Szegedy 2004（walk 算符的谱）；Magniez–Nayak–Roland–Santha 2011
（Search via quantum walk，hitting time 与步数标度）。
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import cast

from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.contracts import (
    OracleView,
    positive_integer,
    require_instance,
    validate_signature,
)
from oracq.algorithms.input_model.interfaces import (
    StatePreparationProtocol,
    checked_state_preparation,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    StatePreparation,
    XorDatabase,
    annotate,
    declare,
    invoke,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import QRAM, Bits, Ref, ValidationError, fuse


@dataclass(frozen=True)
class AdjacencyOracle(OracleView):
    """邻居表查询：``|v>|j>|c> -> |v>|j>|c XOR N(v,j)>``。

    出边表按顶点编号排列，每行恰好 D=2**degree_bits 项；不足处由调用方以自环等方式补齐，
    以保证邻居叠加 1/sqrt(D) 的归一化。paradigm 复用 database_xor，
    address 视为 vertex 与 index 的低位拼接。"""

    oracle_kind = "database_xor"
    operation: Operation

    def adjacency_oracle(self) -> AdjacencyOracle:
        """返回自身的 ``AdjacencyOracle`` 视图，作为邻接 oracle 协议的适配入口。

        Returns:
            AdjacencyOracle: 该视图自身。
        """
        return self

    def __post_init__(self) -> None:
        """校验操作签名以及 vertex/neighbor 等宽等结构约束。"""
        validate_signature(self.operation, ("vertex", "index", "neighbor"), "AdjacencyOracle")
        regs = {r.name: r.type for r in self.operation.module.registers}
        if regs["vertex"].width != regs["neighbor"].width:
            raise ValidationError("vertex 与 neighbor 位宽必须一致")
        if regs["index"].width < 0:
            raise ValidationError("index 位宽无效")

    @property
    def vertex_bits(self) -> int:
        """顶点寄存器 ``vertex`` 的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "vertex")

    @property
    def degree_bits(self) -> int:
        """出边下标寄存器 ``index`` 的位宽。"""
        return next(r.type.width for r in self.operation.module.registers if r.name == "index")

    @property
    def vertices(self) -> int:
        """可寻址的顶点总数，即 ``2**vertex_bits``。"""
        return 1 << self.vertex_bits

    @property
    def degree(self) -> int:
        """补齐后每个顶点的出边表列数 ``D``，即 ``2**degree_bits``。"""
        return 1 << self.degree_bits

    def xor_database(self) -> XorDatabase:
        """把邻居表视为 address=``vertex||index``、data=neighbor 的 XOR database。

        Returns:
            XorDatabase: 包装同一邻居表操作的 XOR 数据库视图。
        """
        from oracq.algorithms.input_model.oracles import XorDatabase

        v, g = self.vertex_bits, self.degree_bits
        b = Builder(
            _name("adjacency_as_database", self.operation),
            {"address": Bits(v + g), "data": Bits(v)},
            resources_for(("adj", self.operation)),
        )
        invoke(
            b,
            self.operation,
            "adj",
            vertex=b["address"][:v],
            index=b["address"][v:],
            neighbor=b["data"],
        )
        return XorDatabase(annotate(b.finish(), "database_xor"))


def as_adjacency(value: AdjacencyOracle | Operation) -> AdjacencyOracle:
    """把图邻接输入适配为 ``AdjacencyOracle`` 视图。

    Args:
        value: 已是 ``AdjacencyOracle`` 时原样返回；否则须是带
            ``(vertex, index, neighbor)`` 签名的完整 ``Operation``。

    Returns:
        AdjacencyOracle: 包装给定操作的邻接 oracle 视图。

    Raises:
        ValidationError: value 既不是 ``AdjacencyOracle`` 也不是 ``Operation``。
    """
    if isinstance(value, AdjacencyOracle):
        return value
    return AdjacencyOracle(require_instance(value, Operation, "adjacency"))


def _as_operation(value: object, path: str) -> Operation:
    """把视图或裸 ``Operation`` 统一为经过程序校验的 ``Operation``。"""
    operation = getattr(value, "operation", None)
    if operation is None:
        operation = value
    require_instance(operation, Operation, path)
    cast("Operation", operation).program()
    return cast("Operation", operation)


def abstract_adjacency(
    vertex_bits: int, degree_bits: int, *, name: str | None = None
) -> AdjacencyOracle:
    """图邻接的开放声明；体为空，可经 bind 绑定 gate/QRAM 实现。

    Args:
        vertex_bits: 顶点寄存器位宽，范围 1..32。
        degree_bits: 出边下标寄存器位宽，范围 0..32。
        name: 声明槽的名称；缺省为 ``GraphAdjacency_{vertex_bits}_{degree_bits}``。

    Returns:
        AdjacencyOracle: 体为空、待绑定实现的邻接声明。
    """
    positive_integer(vertex_bits, "abstract_adjacency.vertex_bits", maximum=32)
    positive_integer(degree_bits, "abstract_adjacency.degree_bits", minimum=0, maximum=32)
    return AdjacencyOracle(
        declare(
            name or f"GraphAdjacency_{vertex_bits}_{degree_bits}",
            {
                "vertex": Bits(vertex_bits),
                "index": Bits(degree_bits),
                "neighbor": Bits(vertex_bits),
            },
            paradigm="database_xor",
            attributes={"oracle_role": "graph_adjacency"},
        )
    )


def _check_table(
    neighbors: Sequence[Sequence[int]],
) -> tuple[tuple[tuple[int, ...], ...], int, int]:
    """校验邻居表为非空矩形且表项是合法顶点编号，返回规整后的行表。"""
    rows = tuple(tuple(row) for row in neighbors)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValidationError("邻居表必须是非空矩形")
    n = len(rows)
    for row in rows:
        for target in row:
            if type(target) is not int or not 0 <= target < n:
                raise ValidationError("邻居表项必须是合法顶点编号")
    return rows, n, len(rows[0])


def _widths(n: int, d: int) -> tuple[int, int]:
    """顶点数 ``n`` 与出边列数 ``d`` 各自所需的编码位宽。"""
    return max(1, (n - 1).bit_length()), (d - 1).bit_length()


def gate_adjacency(
    neighbors: Sequence[Sequence[int]], *, name: str | None = None
) -> AdjacencyOracle:
    """小规模邻居表的门级实现：对每个 (v, j) 分支做受控 X 网络。

    Args:
        neighbors: 非空矩形邻居表，行为顶点 0..V-1、每行 D 项，表项为合法
            顶点编号。
        name: 生成操作的名称；缺省由表内容派生。

    Returns:
        AdjacencyOracle: 受控 X 网络实现的邻接 oracle。
    """
    rows, n, d = _check_table(neighbors)
    v, g = _widths(n, d)
    b = Builder(
        name or _name("gate_adjacency", rows),
        {"vertex": Bits(v), "index": Bits(g), "neighbor": Bits(v)},
    )
    for vertex, row in enumerate(rows):
        for index, target in enumerate(row):
            if target:
                with b.control(fuse(b["vertex"], b["index"]), vertex | (index << v)):
                    for bit in range(v):
                        if (target >> bit) & 1:
                            b.x(b["neighbor"][bit])
    return AdjacencyOracle(
        annotate(b.finish(), "database_xor", implementation="gate_neighbor_table")
    )


def qram_adjacency(
    vertex_bits: int, degree_bits: int, *, name: str | None = None
) -> AdjacencyOracle:
    """邻居表的 QRAM 实现；地址为 ``vertex | (index << vertex_bits)``，字为邻居编号。

    Args:
        vertex_bits: 顶点寄存器位宽，范围 1..32。
        degree_bits: 出边下标寄存器位宽，范围 0..32。
        name: 生成操作的名称；缺省为 ``qram_adjacency_{vertex_bits}_{degree_bits}``。

    Returns:
        AdjacencyOracle: 以 table bank 为后端的邻接 oracle。
    """
    positive_integer(vertex_bits, "qram_adjacency.vertex_bits", maximum=32)
    positive_integer(degree_bits, "qram_adjacency.degree_bits", minimum=0, maximum=32)
    b = Builder(
        name or f"qram_adjacency_{vertex_bits}_{degree_bits}",
        {
            "vertex": Bits(vertex_bits),
            "index": Bits(degree_bits),
            "neighbor": Bits(vertex_bits),
        },
        {"table": QRAM(vertex_bits + degree_bits, vertex_bits)},
    )
    b.qram("table", fuse(b["vertex"], b["index"]), b["neighbor"])
    return AdjacencyOracle(
        annotate(b.finish(), "database_xor", implementation="qram_neighbor_table")
    )


def _neighbor_reflection(
    b: Builder, adjacency: AdjacencyOracle, prefix: str, vertex: Ref, peer: Ref, index: Ref
) -> None:
    """关于 ``span_v{ A|v>|0,0> }`` 的反射，``A|v>|0,0> = 1/sqrt(D) sum_j |v>|N(v,j)>|j>``。

    邻居态制备 A = Adj · H_index 整体上酉（Adj 自逆），故 ``R = A(2|0><0|-I)A†``
    可按 Adj, H, Phi, H, Adj 顺序组装；index 寄存器留在行走空间内。"""
    invoke(b, adjacency.operation, prefix, vertex=vertex, index=index, neighbor=peer)
    b.h(index)
    reflect_zero(b, fuse(peer, index), positive=True)
    b.h(index)
    invoke(b, adjacency.operation, prefix, vertex=vertex, index=index, neighbor=peer)


def szegedy_walk(
    adjacency: AdjacencyOracle | Operation, *, name: str | None = None
) -> Operation:
    """Szegedy 行走步 W = R_B R_A，作用在 (current, peer, index) 行走空间上。

    R_A 是关于"current 端邻居叠加"子空间的反射，R_B 交换两端角色（等价于 SWAP 共轭）。
    图不需要正则；非正则时由补齐方式决定有效转移。

    Args:
        adjacency: 邻接 oracle，或带 ``(vertex, index, neighbor)`` 签名的操作。
        name: 生成操作的名称；缺省由邻接操作派生。

    Returns:
        Operation: 一步行走的酉操作。
    """
    adjacency = as_adjacency(adjacency)
    v, g = adjacency.vertex_bits, adjacency.degree_bits
    b = Builder(
        name or _name("szegedy_walk", adjacency.operation),
        {"current": Bits(v), "peer": Bits(v), "index": Bits(g)},
        resources_for(("adj", adjacency.operation)),
        attributes={
            "algorithm": "szegedy_walk",
            "vertex_bits": v,
            "degree_bits": g,
            "composition": "R_B.R_A",
        },
    )
    _neighbor_reflection(b, adjacency, "adj", b["current"], b["peer"], b["index"])
    _neighbor_reflection(b, adjacency, "adj", b["peer"], b["current"], b["index"])
    return b.finish()


def szegedy_setup(
    adjacency: AdjacencyOracle | Operation, *, name: str | None = None
) -> StatePreparation:
    """制备 MNRS 初态 ``1/sqrt(N) sum_v A|v>``：current 均匀叠加后接邻居态制备。

    target 布局与 szegedy_walk 的寄存器声明顺序一致：低 vertex_bits 位是 current，
    接下来 vertex_bits 位是 peer，高 degree_bits 位是 index。

    Args:
        adjacency: 邻接 oracle，或带 ``(vertex, index, neighbor)`` 签名的操作。
        name: 生成操作的名称；缺省由邻接操作派生。

    Returns:
        StatePreparation: 制备 MNRS 初态的制备视图。
    """
    adjacency = as_adjacency(adjacency)
    v, g = adjacency.vertex_bits, adjacency.degree_bits
    b = Builder(
        name or _name("szegedy_setup", adjacency.operation),
        {"target": Bits(2 * v + g), "work": Bits(0)},
        resources_for(("adj", adjacency.operation)),
    )
    current, peer, index = b["target"][:v], b["target"][v : 2 * v], b["target"][2 * v :]
    b.h(current)
    b.h(index)
    invoke(b, adjacency.operation, "adj", vertex=current, index=index, neighbor=peer)
    return StatePreparation(
        annotate(
            b.finish(),
            "state_prep_isometry",
            zero_input=True,
            clean_work=True,
            implementation="uniform_vertex_plus_neighbor",
        )
    )


def quantum_walk_search(
    setup: StatePreparationProtocol,
    walk: Operation,
    marked: Operation,
    steps: int,
    *,
    name: str | None = None,
) -> Operation:
    """MNRS marked-vertex 搜索骨架：setup 后 Repeat{ walk; marked 相位翻转 }。

    Args:
        setup: 初态制备（StatePreparation 协议），target 为整个行走空间。
        walk: 行走步句柄；其寄存器按声明顺序映射到 setup.target 的连续切片。
        marked: 相位 oracle（target[, work] 接口），作用于行走空间的前 target.width 位，对 szegedy_walk 而言即 current 寄存器。
        steps: 非负迭代次数，可参照 suggest_steps 选择。

    三个句柄都允许是开放声明；生成的程序经 bind 分批绑定实现后即可执行。

    Returns:
        Operation: 重复行走与相位翻转的搜索程序。
    """
    positive_integer(steps, "quantum_walk_search.steps", minimum=0)
    prep = checked_state_preparation(setup, path="quantum_walk_search.setup")
    walk_op = _as_operation(walk, "quantum_walk_search.walk")
    marked_op = _as_operation(marked, "quantum_walk_search.marked")
    width = sum(r.type.width for r in walk_op.module.registers)
    if width != prep.width:
        raise ValidationError("walk 与 setup 的行走空间宽度不一致")
    marked_regs = {r.name: r.type for r in marked_op.module.registers}
    if "target" not in marked_regs or set(marked_regs) - {"target", "work"}:
        raise ValidationError("marked 相位 oracle 需要 target[, work] 接口")
    marked_width = marked_regs["target"].width
    marked_work = marked_regs.get("work", Bits(0)).width
    if marked_width > width:
        raise ValidationError("marked 目标超出行走空间宽度")
    b = Builder(
        name or _name("quantum_walk_search", prep.operation, walk_op, marked_op, steps),
        {"target": Bits(width), "work": Bits(prep.work_width + marked_work)},
        resources_for(("setup", prep.operation), ("walk", walk_op), ("marked", marked_op)),
        attributes={
            "algorithm": "quantum_walk_search",
            "framework": "MNRS",
            "steps": steps,
            "walk_width": width,
        },
    )
    setup_work, marked_ref = b["work"][: prep.work_width], b["work"][prep.work_width :]
    invoke(b, prep.operation, "setup", target=b["target"], work=setup_work)
    with b.repeat(steps):
        cursor = 0
        arguments: dict[str, Ref] = {}
        for reg in walk_op.module.registers:
            arguments[reg.name] = (
                b["target"][cursor : cursor + reg.type.width].reinterpret(reg.type.kind)
            )
            cursor += reg.type.width
        invoke(b, walk_op, "walk", **arguments)
        marked_arguments = {"target": b["target"][:marked_width]}
        if marked_work:
            marked_arguments["work"] = marked_ref
        invoke(b, marked_op, "marked", **marked_arguments)
    return b.finish()


def transition_matrix(neighbors: Sequence[Sequence[int]]) -> tuple[tuple[float, ...], ...]:
    """由邻居表构造经典随机游走转移矩阵，``P[v][u] = |{j: N(v,j)=u}| / D``。

    Args:
        neighbors: 非空矩形邻居表，表项为合法顶点编号。

    Returns:
        tuple[tuple[float, ...], ...]: 行主序的转移概率矩阵，每行和为 1。
    """
    rows, n, d = _check_table(neighbors)
    result = []
    for row in rows:
        counts = [0.0] * n
        for target in row:
            counts[target] += 1.0 / d
        result.append(tuple(counts))
    return tuple(result)


def hitting_times(
    transition: Sequence[Sequence[float]], marked: Iterable[int]
) -> tuple[float, ...]:
    """每个顶点首次到达 marked 集合的期望步数（marked 顶点为 0）。

    解线性方程 (I - P_free) h = 1，P_free 为删去 marked 行/列的子转移矩阵。

    Args:
        transition: 行和为 1 的非负转移方阵。
        marked: 标记顶点编号集合，须非空且不越界。

    Returns:
        tuple[float, ...]: 按顶点编号排列的期望首达步数。
    """
    n = len(transition)
    if not n or any(len(row) != n for row in transition):
        raise ValidationError("转移矩阵必须是非空方阵")
    for row in transition:
        if any(type(p) not in (int, float) or p < 0 for p in row):
            raise ValidationError("转移概率必须是非负实数")
        if abs(sum(row) - 1.0) > 1e-9:
            raise ValidationError("转移矩阵行和必须为 1")
    marked_set = set(marked)
    if not marked_set or any(type(m) is not int or not 0 <= m < n for m in marked_set):
        raise ValidationError("marked 必须是合法顶点编号的非空集合")
    free = [v for v in range(n) if v not in marked_set]
    if not free:
        return tuple(0.0 for _ in range(n))
    size = len(free)
    system = [
        [
            (1.0 if i == j else 0.0) - float(transition[free[i]][free[j]])
            for j in range(size)
        ]
        + [1.0]
        for i in range(size)
    ]
    for col in range(size):
        pivot = max(range(col, size), key=lambda r: abs(system[r][col]))
        if abs(system[pivot][col]) < 1e-12:
            raise ValidationError("marked 集合从某些顶点不可达，首达时间发散")
        system[col], system[pivot] = system[pivot], system[col]
        for row in range(size):  # type: ignore[assignment]
            if row != col:
                factor = system[cast("int", row)][col] / system[col][col]
                for k in range(col, size + 1):
                    system[cast("int", row)][k] -= factor * system[col][k]
    hits = {free[i]: system[i][size] / system[i][i] for i in range(size)}
    return tuple(0.0 if v in marked_set else hits[v] for v in range(n))


def suggest_steps(transition: Sequence[Sequence[float]], marked: Iterable[int]) -> int:
    """MNRS 建议行走步数：ceil(pi/4 * sqrt(H_avg))，H_avg 为平均首达时间。

    与 Grover 极限一致：完全图上 ``H_avg = N/|M|``，步数恢复 ``pi/4 * sqrt(N/|M|)``。

    Args:
        transition: 行和为 1 的非负转移方阵。
        marked: 标记顶点编号集合，须非空且不越界。

    Returns:
        int: 建议迭代次数，至少为 1；无自由顶点时为 0。
    """
    hits = hitting_times(transition, marked)
    n = len(transition)
    free = n - len(set(marked))
    if not free:
        return 0
    average = sum(hits) / free
    return max(1, math.ceil(math.pi / 4 * math.sqrt(average)))
