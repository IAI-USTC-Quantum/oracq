"""QFVM 数据路径的 QMem 直连重写：状态表/几何表经指针式二维寻址，残差树用 QVector。

电路结构与 applications/qfvm.py 同构（可逆 Roe 算术、九槽位几何、补齐对角、
原地位置置换全部保持），但不再经过抽象数据库槽位与 bind：模块直接声明
QRAM 形式资源，三守恒量合并为一张 (场, 单元) 状态表，邻居 cell±1 用模 2^cell_width
的指针算术实现周期边界，几何表按 (列, 槽位) 二维寻址，残差态由 QVector
（平方范数树 + 符号反冲）制备。语义与既有路径等价性由 tests/core/test_qfvm_qmem.py
在 simulate 上逐振幅对拍。
"""

from __future__ import annotations

from pyqecclang.algorithms.common.arithmetic import BooleanNetwork, FixedFormat
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import annotate
from pyqecclang.algorithms.input_model.qdata import QVector
from pyqecclang.algorithms.input_model.sparse import compare_words, value_transposition
from pyqecclang.applications.flow_data import RoeFlowData
from pyqecclang.applications.qfvm import geometry_cells
from pyqecclang.applications.roe import ArithmeticBuilder, roe_face
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import QRAM, Adjoint, Bits, ValidationError, fuse
from pyqecclang.infrastructure.qmem import QMem


def _geometry_refs(ref, inputs):
    sizes = (inputs.width, 4, inputs.cell_width, 2, 2, 2, 1)
    refs, offset = [], 0
    for size in sizes:
        refs.append(ref[offset : offset + size])
        offset += size
    return refs


def qfvm_qmem_physical(inputs, *, gamma=1.4, entropy_delta=0.125, mass=1.0, dx=1.0):
    """矩阵元物理计算：状态表按 (场, 邻居单元) 二维寻址，周期邻居用模加指针。"""
    cw, fmt = inputs.cell_width, inputs.fmt
    b = Builder(
        _name("qfvm_qmem_physical", cw, fmt, gamma, entropy_delta, mass, dx),
        {
            "source": Bits(cw),
            "row": Bits(2),
            "col": Bits(2),
            "band": Bits(2),
            "value": Bits(fmt.width),
            "status": Bits(2),
        },
        {"state": QRAM(cw + 2, fmt.width)},
        attributes={
            "algorithm": "qfvm_arithmetic_entry",
            "correctness": "pending",
            "data_access": "qmem_state_table",
            "boundary": "periodic",
            "mass": mass,
            "dx": dx,
        },
    )
    state = QMem(b, "state", shape=(3, 1 << cw))
    g = ArithmeticBuilder(b, fmt)
    words = []
    for offset in (-1, 0, 1):
        addr = g.local(cw)
        b.xor(b["source"], addr)
        b.add_const(addr.reinterpret("uint"), offset % (1 << cw))
        row_words = []
        for field in range(3):
            word = g.local()
            state[field, addr].load(word)
            row_words.append(word)
        words.append(row_words)
    face = roe_face(fmt=fmt, gamma=gamma, entropy_delta=entropy_delta)
    faces = []
    for i in range(2):
        left, right, flag = g.local(), g.local(), g.local(2)
        b.call(
            face,
            **dict(
                zip(
                    ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r"),
                    words[i] + words[i + 1],
                    strict=True,
                )
            ),
            row=b["row"],
            col=b["col"],
            left=left,
            right=right,
            status=flag,
        )
        g.flags.append(flag)
        faces.append((g.word(left), g.word(right)))
    equal = b.local("component_equal", Bits(1))
    net = BooleanNetwork()
    r, c = net.input("r", 2), net.input("c", 2)
    net.outputs = {"equal": [net.inv(net.any([net.xor(x, y) for x, y in zip(r, c, strict=True)]))]}
    b.call(net.operation(), r=b["row"], c=b["col"], equal=equal)
    west = -faces[0][0] / dx
    center = (faces[1][0] - faces[0][1]) / dx + g.choose(equal, mass, 0)
    east = faces[1][1] / dx
    value = g.choose3(b["band"], [west, center, east])
    return g.finish([(b["value"], value)], b["status"])


def qfvm_qmem_location(inputs):
    """CKS 原地位置置换：九个结构槽位经 (列, 槽位) 二维几何寻址。"""
    n = inputs.width
    b = Builder(
        _name("qfvm_qmem_location", n),
        {"column": Bits(n), "index": Bits(n), "work": Bits(n)},
        {"geometry": QRAM(n + 4, inputs.geometry_width)},
        attributes={
            "sparsity": 9,
            "orientation": "column",
            "padding": "identity on variable 3",
            "full_permutation_extension": True,
            "construction": "nine coherent transpositions",
            "data_access": "qmem_geometry_2d",
        },
    )
    geometry = QMem(b, "geometry", shape=(16, 1 << n))
    neighbors, lefts = [], []
    for rank in range(9):
        slot = b.local("slot_" + str(rank), Bits(4))
        data = b.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        neighbor = b.local("neighbor_" + str(rank), Bits(n))
        left = b.local("left_" + str(rank), Bits(n))
        for bit in range(4):
            if (rank >> bit) & 1:
                b.x(slot[bit])
        geometry[slot, b["column"]].load(data)
        b.xor(data[:n], neighbor)
        padded = b.local("padded_neighbor_" + str(rank), Bits(n))
        b.xor(b["column"], padded)
        b.add_const(padded.reinterpret("uint"), rank)
        with b.control(b["column"][:2], 3):
            b.xor(data[:n], neighbor)
            b.xor(padded, neighbor)
        for bit in range(n):
            if (rank >> bit) & 1:
                b.x(left[bit])
        for previous in range(rank):
            b.call(value_transposition(n), index=left, a=lefts[previous], b=neighbors[previous])
        neighbors.append(neighbor)
        lefts.append(left)
    setup = tuple(b._frames[0])
    for left, right in zip(lefts, neighbors, strict=True):
        b.call(value_transposition(n), index=b["index"], a=left, b=right)
    b.emit(Adjoint(setup))
    return annotate(b.finish(), "sparse_location_inplace")


def qfvm_qmem_entry(inputs, *, padding_value=1.0, **entry_options):
    """任意坐标矩阵元 oracle：与 qfvm_sparse_access 的 entry 同语义，数据面直连 QMem。"""
    if padding_value <= 0 or inputs.fmt.decode(inputs.fmt.encode(padding_value)) != padding_value:
        raise ValidationError("补齐对角值必须是定点格式可精确表示的正数")
    n = inputs.width
    physical = qfvm_qmem_physical(inputs, **entry_options)
    b = Builder(
        _name("qfvm_qmem_entry", n, inputs.fmt, padding_value),
        {"row": Bits(n), "column": Bits(n), "data": Bits(inputs.fmt.width)},
        {
            "geometry": QRAM(n + 4, inputs.geometry_width),
            "state": QRAM(inputs.cell_width + 2, inputs.fmt.width),
        },
        attributes={
            "value_encoding": "signed_fixed_point",
            "value_fraction": inputs.fmt.fraction,
            "matrix": "Hermitian dilation of Roe matrix, identity on padded coordinates",
            "invalid_arithmetic": "entry totalized to zero",
            "correctness": "pending",
            "data_access": "qmem_state_geometry",
        },
    )
    geometry = QMem(b, "geometry", shape=(16, 1 << n))
    selected = b.local("selected_geometry", Bits(inputs.geometry_width))
    for rank in range(9):
        slot = b.local("slot_" + str(rank), Bits(4))
        geom = b.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        match = b.local("match_" + str(rank), Bits(1))
        for bit in range(4):
            if (rank >> bit) & 1:
                b.x(slot[bit])
        geometry[slot, b["column"]].load(geom)
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        with b.control(fuse(match, geom[inputs.geometry_width - 1])):
            b.xor(geom[: inputs.geometry_width - 1], selected[: inputs.geometry_width - 1])
            b.x(selected[inputs.geometry_width - 1])
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        geometry[slot, b["column"]].load(geom)
    _, _, source, row, col, band, valid = _geometry_refs(selected, inputs)
    value = b.local("computed_value", Bits(inputs.fmt.width))
    status = b.local("arithmetic_status", Bits(2))
    b.call(
        physical,
        source=source,
        row=row,
        col=col,
        band=band,
        value=value,
        status=status,
        resources={"state": "state"},
    )
    same = b.local("diagonal", Bits(1))
    b.call(compare_words(n), a=b["row"], b=b["column"], flag=same)
    forward = tuple(b._frames[0])
    with b.control(valid):
        with b.control(status, 0):
            b.xor(value, b["data"])
    with b.control(fuse(same, b["column"][:2]), 7):
        raw = inputs.fmt.encode(padding_value)
        for bit in range(inputs.fmt.width):
            if (raw >> bit) & 1:
                b.x(b["data"][bit])
    b.emit(Adjoint(forward))
    return annotate(b.finish(), "sparse_entry_xor")


class RoeQmemData:
    """经典侧流场数据：复用 RoeFlowData 的物理量与残差，量子面重组为 QMem bank。"""

    def __init__(self, states, *, fmt=None, angle_width=10, gamma=1.4, entropy_delta=0.125, dx=1.0,
                 name="residual"):
        fmt = fmt or FixedFormat(10, 5)
        self.flow = RoeFlowData(
            states, fmt=fmt, angle_width=angle_width, gamma=gamma,
            entropy_delta=entropy_delta, dx=dx,
        )
        self.fmt, self.angle_width = fmt, angle_width
        self.cell_width = self.flow.n.bit_length() - 1
        values = []
        for cell in range(self.flow.n):
            residual = self.flow.residuals[cell]
            values.extend(residual[j] if j < 3 else 0.0 for j in range(4))
        self.vector = QVector(values, fmt=fmt, angle_width=angle_width, name=name)

    @property
    def state_bank(self):
        """(场, 单元) 状态表：地址 = field·n + cell，字为守恒量定点编码。"""
        banks = self.flow.store.snapshot()
        table = {}
        for field, key in enumerate(("rho", "momentum", "energy")):
            for address, word in banks[key].items():
                table[(field << self.cell_width) | address] = word
        return table

    def memories(self, inputs):
        """导出 QMem 直连数据路径的全部运行时内存表。

        Args:
            inputs: ``roe_qfvm_inputs`` 声明的槽位集合，决定几何表布局。

        Returns:
            dict: ``state`` 状态表、``geometry`` 几何表，以及残差 ``QVector``
            快照的角字与符号 bank（默认 ``residual_angles`` 和
            ``residual_sign``，名字随构造时的 ``name`` 参数）；供执行入口
            按名绑定模块声明的 QRAM 资源。
        """
        return {
            "state": self.state_bank,
            "geometry": geometry_cells(inputs),
            **self.vector.snapshot(signed=True),
        }


def qfvm_qmem_rhs(data: RoeQmemData):
    """残差态制备：QVector 平方范数树 + 符号相位反冲（替代 qram_state_prep+sign 组合）。"""
    return data.vector.preparation(signed=True)
