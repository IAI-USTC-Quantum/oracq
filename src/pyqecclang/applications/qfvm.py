"""从 QRAM 流场与 Roe 算术构造 QFVM 稀疏输入及可替换 QLSS 问题。"""

from __future__ import annotations

import math
from dataclasses import dataclass

from pyqecclang.algorithms.arithmetic import FixedFormat
from pyqecclang.algorithms.operators import _name
from pyqecclang.algorithms.oracles import (
    SparseAccess,
    StatePreparation,
    abstract_database,
    abstract_state_prep,
    annotate,
    invoke,
    qram_database,
    qram_state_prep,
    resources_for,
)
from pyqecclang.algorithms.sparse import compare_words, value_transposition
from pyqecclang.applications.roe import ArithmeticBuilder, roe_face
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Adjoint, Bits, ValidationError, fuse
from pyqecclang.infrastructure.linking import Binding, bind


@dataclass(frozen=True)
class RoeQfvmInputs:
    """QFVM 输入模型声明的抽象槽位集合：原始场量、几何、残差值与残差态制备。

    槽位本身不含数据；实现由 ``bind_qfvm`` 绑定，运行时内存表由 ``qfvm_memories``
    物化。矩阵元素不在任何槽位预存，而由可逆 Roe 算术在查询时现场计算。

    Attributes:
        cell_width: 单元编号位宽，网格单元数为 ``2**cell_width``。
        fmt: 守恒量与残差值使用的定点格式。
        angle_width: 旋转角字的位宽。
        rho: 密度（第一个守恒量）数据库的抽象声明槽位。
        momentum: 动量（第二个守恒量）数据库的抽象声明槽位。
        energy: 能量（第三个守恒量）数据库的抽象声明槽位。
        geometry: 几何表数据库的槽位，只编码稀疏位置与索引，不含矩阵值。
        theta: 残差值定点字到旋转角字换算表的槽位（预留数据面）。
        rhs: 归一化残差态制备的抽象槽位。
        residual: 每 (单元, 分量) 量化残差值数据库的槽位。
    """

    cell_width: int
    fmt: FixedFormat
    angle_width: int
    rho: object
    momentum: object
    energy: object
    geometry: object
    theta: object
    rhs: StatePreparation
    residual: object

    @property
    def width(self):
        """完整矩阵坐标位宽：2 位分量 + cell_width 位单元 + 1 位扩张半块标志。"""
        return self.cell_width + 3

    @property
    def geometry_width(self):
        """几何表字位宽：对手坐标、对称槽位、源单元、行/列分量、源带与有效位七段打包。"""
        return self.width + 4 + self.cell_width + 2 + 2 + 2 + 1


def roe_qfvm_inputs(*, cell_width=2, fmt=None, angle_width=8, prefix="RoeQfvm"):
    """声明 QFVM 的整套抽象输入槽位：六个数据库加一个 RHS 态制备。

    返回的全是开放槽位，尚无数据；数据到绑定（``bind_qfvm``）与执行期
    （``qfvm_memories``）才出现。

    Args:
        cell_width: 单元编号位宽，范围为 2..25。
        fmt: 守恒量定点格式，省略时使用 ``FixedFormat(6, 2)``。
        angle_width: 旋转角字位宽，范围为 1..64。
        prefix: 各抽象槽位模块名的前缀。

    Returns:
        RoeQfvmInputs: 声明集合，含三守恒量库、几何表、theta 表、残差值库
        与抽象 RHS 制备。

    Raises:
        ValidationError: cell_width 或 angle_width 超出范围。
    """
    fmt = fmt or FixedFormat(6, 2)
    if not 2 <= cell_width <= 25 or not 1 <= angle_width <= 64:
        raise ValidationError("QFVM 需要至少四个单元，且几何/角度字长不超过 64 位")
    width = cell_width + 3
    geometry_width = width + 4 + cell_width + 7
    return RoeQfvmInputs(
        cell_width,
        fmt,
        angle_width,
        abstract_database(prefix + "Rho", cell_width, fmt.width),
        abstract_database(prefix + "Momentum", cell_width, fmt.width),
        abstract_database(prefix + "Energy", cell_width, fmt.width),
        abstract_database(prefix + "Geometry", width + 4, geometry_width),
        abstract_database(prefix + "Theta", fmt.width, angle_width),
        abstract_state_prep(prefix + "Rhs", width, cell_width + 2 + angle_width + 1),
        abstract_database(prefix + "RhsValues", cell_width + 2, fmt.width),
    )


def geometry_cells(inputs):
    """只编码稀疏位置与原始数据索引，不含任何流场矩阵值。"""
    cw, w = inputs.cell_width, inputs.width
    n, table = 1 << cw, {}
    for row in range(1 << w):
        cell, var, half = (row >> 2) % n, row % 4, row >> (cw + 2)
        for slot in range(16):
            band, other_var = divmod(slot, 3)
            valid = int(slot < 9 and var < 3)
            neighbor_cell = (cell + band - 1) % n
            neighbor = other_var + 4 * neighbor_cell + ((1 - half) << (cw + 2))
            reverse = (2 - band) * 3 + var if valid else 0
            source = cell if half == 0 else neighbor_cell
            rowvar, colvar, source_band = (
                (var, other_var, band) if half == 0 else (other_var, var, 2 - band)
            )
            values = (
                (neighbor, w),
                (reverse, 4),
                (source, cw),
                (rowvar, 2),
                (colvar, 2),
                (source_band, 2),
                (valid, 1),
            )
            packed, offset = 0, 0
            for value, width in values:
                packed |= (value & ((1 << width) - 1)) << offset
                offset += width
            table[row + (slot << w)] = packed
    return table


def ptheta_cells(fmt, angle_width, amax):
    """生成 theta 表：残差值定点字到旋转角字 ``2*acos(min(1, |v|/amax))`` 的换算。

    Args:
        fmt: 定点格式，地址取该格式的全部 ``2**fmt.width`` 个原始字。
        angle_width: 输出角字位宽。
        amax: 元素幅值上界，用于截断 acos 的自变量。

    Returns:
        dict: 原始定点字到量化角字的映射，角字按 ``2**angle_width`` 回绕。

    Raises:
        ValueError: amax 非正。
    """
    if amax <= 0:
        raise ValueError("Amax 必须为正")
    return {
        raw: round(
            2 * math.acos(min(1, abs(fmt.decode(raw)) / amax)) * (1 << angle_width) / (2 * math.pi)
        )
        % (1 << angle_width)
        for raw in range(1 << fmt.width)
    }


def roe_entry(inputs, *, gamma=1.4, entropy_delta=0.125, mass=1.0, dx=1.0):
    """构造单个 Roe 矩阵元的可逆算术电路（不含稀疏寻址）。

    对 ``source`` 及其左右邻居单元各查询三个守恒量库（周期边界回绕），
    调用 ``roe_face`` 计算两个界面的通量特征分量，按 ``band`` 组装西/中/东
    三个候选值；质量项只在行列分量相同时加到中心块对角。

    Args:
        inputs: ``roe_qfvm_inputs`` 声明的槽位集合，提供守恒量库与定点格式。
        gamma: 比热比。
        entropy_delta: Roe 熵修正系数。
        mass: 质量项，只作用于中心块的分量对角。
        dx: 网格步长。

    Returns:
        Operation: 寄存器为 source、row、col、band、value、status；value 写出
        选定矩阵元，status 聚合各算术节点的失效旗标（非零表示溢出或除零）。
    """
    cw, fmt = inputs.cell_width, inputs.fmt
    fields = (
        ("rho", inputs.rho.operation),
        ("momentum", inputs.momentum.operation),
        ("energy", inputs.energy.operation),
    )
    b = Builder(
        _name("roe_entry", cw, fmt, gamma, entropy_delta, mass, dx, *[o for _, o in fields]),
        {
            "source": Bits(cw),
            "row": Bits(2),
            "col": Bits(2),
            "band": Bits(2),
            "value": Bits(fmt.width),
            "status": Bits(2),
        },
        resources_for(*fields),
        attributes={
            "algorithm": "qfvm_arithmetic_entry",
            "correctness": "pending",
            "matrix_table": False,
            "boundary": "periodic",
            "mass": mass,
            "dx": dx,
        },
    )
    g = ArithmeticBuilder(b, fmt)
    addresses, words = [], []
    for offset in (-1, 0, 1):
        addr = g.local(cw)
        b.xor(b["source"], addr)
        b.add_const(addr.reinterpret("uint"), offset % (1 << cw))
        addresses.append(addr)
        state = []
        for name, op in fields:
            word = g.local()
            invoke(b, op, name, address=addr, data=word)
            state.append(word)
        words.append(state)
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
    # 质量项只作用于中心块的分量对角。
    equal = b.local("component_equal", Bits(1))
    from pyqecclang.algorithms.arithmetic import BooleanNetwork

    net = BooleanNetwork()
    r, c = net.input("r", 2), net.input("c", 2)
    net.outputs = {"equal": [net.inv(net.any([net.xor(x, y) for x, y in zip(r, c, strict=True)]))]}
    b.call(net.operation(), r=b["row"], c=b["col"], equal=equal)
    west = -faces[0][0] / dx
    center = (faces[1][0] - faces[0][1]) / dx + g.choose(equal, mass, 0)
    east = faces[1][1] / dx
    value = g.choose3(b["band"], [west, center, east])
    return g.finish([(b["value"], value)], b["status"])


def _geometry_refs(ref, inputs):
    sizes = (inputs.width, 4, inputs.cell_width, 2, 2, 2, 1)
    refs, offset = [], 0
    for size in sizes:
        refs.append(ref[offset : offset + size])
        offset += size
    return refs


def roe_qfvm_block_encoding(inputs, *, amax=8.0, padding_value=1.0, **entry_options):
    """显式从稀疏输入转换到 BE；QFVM 本身不再强制只暴露 BE。"""
    from pyqecclang.algorithms.sparse import real_symmetric_sparse_encoding

    if not 0 < padding_value <= amax:
        raise ValidationError("补齐对角值必须为正且不超过元素上界")
    access = qfvm_sparse_access(inputs, padding_value=padding_value, **entry_options)
    return real_symmetric_sparse_encoding(access, inputs.fmt, amax, diagonal_nonnegative=True)


def rhs_qram_preparation(inputs):
    """构造归一化残差态的 QRAM 制备：幅度取自平方范数树，符号经相位反冲写入。

    幅度由 ``qram_state_prep`` 按层查询 ``rhs_angles`` 角字 bank 得到；
    符号用一位 ``rhs_sign`` bank 经 Load、Z 与反 Load 把负残差分量变成
    pi 相位。输入为零态，work 复净。

    Returns:
        StatePreparation: target 为完整矩阵坐标（cell_width+3 位），幅度
        只写在前 cell_width+2 位上，扩张半块标志保持零。
    """
    n, aw = inputs.cell_width + 2, inputs.angle_width
    prep = qram_state_prep(n, aw)
    sign = qram_database(n, 1)
    b = Builder(
        _name("roe_rhs", n, aw),
        {"target": Bits(n + 1), "work": Bits(n + aw + 1)},
        resources_for(("prep", prep.operation), ("sign", sign.operation)),
        attributes={"correctness": "pending", "algorithm": "qram_signed_residual_tree"},
    )
    invoke(b, prep.operation, "prep", target=b["target"][:n], work=b["work"][: n + aw])
    flag = b["work"][n + aw :]
    invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)
    b.z(flag)
    invoke(b, sign.operation, "sign", address=b["target"][:n], data=flag)
    return StatePreparation(annotate(b.finish(), "state_prep_isometry", zero_input=True))


def qram_bindings(inputs):
    """给出 QFVM 抽象槽位到 QRAM 实现的完整绑定映射。

    三个守恒量库、几何表与 theta 表绑定 ``qram_database`` 实际库（table
    资源映射到同名经典 bank）；残差值槽绑定 ``RoeResidualValues`` 库
    （rhs_values）；RHS 制备槽绑定 ``rhs_qram_preparation`` 电路（角字与
    符号资源分别映射到 rhs_angles 与 rhs_sign）。

    Returns:
        dict: 抽象槽位模块名到 ``Binding`` 的映射，可直接交给 ``bind``。
    """
    mapping = {}
    for key in ("rho", "momentum", "energy", "geometry", "theta"):
        db = getattr(inputs, key)
        actual = qram_database(db.address_width, db.data_width, name="RoeMemory_" + key)
        mapping[db.operation.module.name] = Binding(actual.operation, {"table": key})
    residual = qram_database(
        inputs.residual.address_width, inputs.residual.data_width, name="RoeResidualValues"
    )
    mapping[inputs.residual.operation.module.name] = Binding(
        residual.operation, {"table": "rhs_values"}
    )
    prep = rhs_qram_preparation(inputs)
    mapping[inputs.rhs.operation.module.name] = Binding(
        prep.operation, {"prep__angles": "rhs_angles", "sign__table": "rhs_sign"}
    )
    return mapping


def bind_qfvm(program, inputs):
    """把程序中尚未解析的 QFVM 抽象槽位绑定到 QRAM 实现与符号残差树制备。

    Args:
        program: 含 QFVM 抽象槽调用的程序。
        inputs: ``roe_qfvm_inputs`` 声明的槽位集合。

    Returns:
        Program: 只替换未解析槽位后的程序；与 QFVM 无关的开放槽保持原样。
    """
    from pyqecclang.infrastructure.linking import unresolved

    missing = {r.name for r in unresolved(program)}
    return bind(program, {k: v for k, v in qram_bindings(inputs).items() if k in missing})


def qfvm_memories(inputs, flow, *, amax=8.0):
    """物化 QFVM 路径的运行时内存表：流场快照加静态几何表与 theta 表。

    Args:
        inputs: ``roe_qfvm_inputs`` 声明的槽位集合。
        flow: 经典侧流场数据，提供六个 bank 的快照。
        amax: theta 表使用的元素幅值上界。

    Returns:
        dict: 按名交给执行后端的内存表：rho、momentum、energy、rhs_values、
        rhs_sign、rhs_angles、geometry 与 theta。

    Raises:
        ValidationError: 流场与声明的单元数、定点或角度格式不匹配。
    """
    if (
        flow.n != (1 << inputs.cell_width)
        or flow.fmt != inputs.fmt
        or flow.angle_width != inputs.angle_width
    ):
        raise ValidationError("QFVM 流场和 oracle 的单元数/定点/角度格式不匹配")
    result = flow.store.snapshot()
    result["geometry"] = geometry_cells(inputs)
    result["theta"] = ptheta_cells(inputs.fmt, inputs.angle_width, amax)
    return result


def roe_qfvm_problem(
    inputs, *, spectrum, rhs_norm=None, amax=8.0, padding_value=1.0, **entry_options
):
    """D=([[0,M],[M.T,0]] on physical coordinates) + padding_value*I_pad。"""
    from pyqecclang.algorithms.qlss import LinearSystem, SparseSystem, SpectralPromise

    if not 0 < padding_value <= amax:
        raise ValidationError("补齐对角值必须为正且不超过元素上界")
    if not isinstance(spectrum, SpectralPromise):
        raise ValidationError("QFVM QLSS 输入需要显式 SpectralPromise")
    extended = SpectralPromise(
        max(spectrum.norm_upper, padding_value),
        min(spectrum.sigma_min_lower, padding_value),
        spectrum.evidence,
    )
    sparse = SparseSystem(
        qfvm_sparse_access(inputs, padding_value=padding_value, **entry_options),
        inputs.fmt,
        amax,
        inputs.rhs,
        extended,
        diagonal_nonnegative=True,
        hermitian=True,
    )
    return LinearSystem(
        sparse=sparse,
        physical_width=inputs.cell_width + 2,
        physical_high_value=1,
        rhs_norm=rhs_norm,
        data_assumptions=(
            "coherent XOR QRAM on arbitrary address superpositions",
            "memory snapshot fixed across coherent calls, adjoints and norm probes",
            "raw conserved fields and O(Ns) geometry; no matrix-entry table",
            "signed residual value oracle and norm tree with cached rotation angles",
            "clean reversible RHS preparation; norm known separately",
            "logical local patches; native QRAM banks currently rematerialized",
        ),
    )


def roe_qfvm_step(inputs, qlss, *, spectrum=None, rhs_norm=None, **options):
    """把 QFVM 问题交给可替换的 QLSS 协议求解。

    Args:
        inputs: ``roe_qfvm_inputs`` 声明的槽位集合。
        qlss: 声明了 ``input_model`` 的 ``QLSSProtocol``。
        spectrum: 覆盖实际量化后矩阵的 ``SpectralPromise``，必填。
        rhs_norm: 经典右端范数；提供后才能经 ``SolveResult.recover_norm``
            恢复解的幅值。
        **options: 其余关键字参数透传给 ``roe_qfvm_problem``，含 amax、
            padding_value 及矩阵元算术选项。

    Returns:
        SolveResult: 协议生成的解态、范数探针与适配记录。

    Raises:
        ValidationError: qlss 不是声明 ``input_model`` 的 ``QLSSProtocol``，
            或未提供 spectrum。
    """
    from pyqecclang.algorithms.qlss import QLSSProtocol

    if not isinstance(qlss, QLSSProtocol):
        raise ValidationError(
            "请使用声明 input_model 的 QLSSProtocol；裸 callable 不能表达输入模型"
        )
    if spectrum is None:
        raise ValidationError(
            "QFVM 替换 QLSS 需要最小奇异值/范数声明：spectrum=SpectralPromise(...)"
        )
    return qlss.solve(roe_qfvm_problem(inputs, spectrum=spectrum, rhs_norm=rhs_norm, **options))


def qfvm_sparse_access(inputs, *, padding_value=1.0, **entry_options):
    """构造 QFVM 的 CKS 稀疏访问：原地位置置换与任意坐标矩阵元 XOR 两个 oracle。

    位置 oracle 查询每列的九个结构槽位，经九次相干值转置补全为完整置换；
    补齐分量没有邻居结构，对手坐标改写为 ``column + rank`` 得到 padding
    对角。条目 oracle 用几何表选择源单元与分量，调用 ``roe_entry`` 现场
    计算元素：结构域外条目为零，算术失效时条目归零，补齐坐标的对角写
    padding_value。

    Args:
        inputs: ``roe_qfvm_inputs`` 声明的槽位集合。
        padding_value: 补齐对角值，须为定点格式可精确表示的正数。
        **entry_options: 透传给 ``roe_entry`` 的算术选项（gamma、mass 等）。

    Returns:
        SparseAccess: 坐标宽 ``inputs.width``、值宽 ``inputs.fmt.width``、
        稀疏度 9。

    Raises:
        ValidationError: padding_value 非正或定点格式无法精确表示。
    """

    if padding_value <= 0 or inputs.fmt.decode(inputs.fmt.encode(padding_value)) != padding_value:
        raise ValidationError("补齐对角值必须是定点格式可精确表示的正数")
    n = inputs.width
    geometry = inputs.geometry.operation
    locator = Builder(
        _name("qfvm_sparse_location", geometry, padding_value),
        {"column": Bits(n), "index": Bits(n), "work": Bits(0)},
        resources_for(("geometry", geometry)),
        attributes={
            "sparsity": 9,
            "orientation": "column",
            "padding": "identity on variable 3",
            "full_permutation_extension": True,
            "construction": "nine coherent transpositions",
        },
    )
    neighbors, lefts = [], []
    for rank in range(9):
        slot = locator.local("slot_" + str(rank), Bits(4))
        data = locator.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        neighbor = locator.local("neighbor_" + str(rank), Bits(n))
        left = locator.local("left_" + str(rank), Bits(n))
        for bit in range(4):
            if (rank >> bit) & 1:
                locator.x(slot[bit])
        invoke(locator, geometry, "geometry", address=fuse(locator["column"], slot), data=data)
        locator.xor(data[:n], neighbor)
        # 补齐变量具有对角 padding_value；其余八个不同位置返回零元素。
        padded = locator.local("padded_neighbor_" + str(rank), Bits(n))
        locator.xor(locator["column"], padded)
        locator.add_const(padded.reinterpret("uint"), rank)
        with locator.control(locator["column"][:2], 3):
            locator.xor(data[:n], neighbor)
            locator.xor(padded, neighbor)
        for bit in range(n):
            if (rank >> bit) & 1:
                locator.x(left[bit])
        for previous in range(rank):
            locator.call(
                value_transposition(n), index=left, a=lefts[previous], b=neighbors[previous]
            )
        neighbors.append(neighbor)
        lefts.append(left)
    setup = tuple(locator._frames[0])
    for left, right in zip(lefts, neighbors, strict=True):
        locator.call(value_transposition(n), index=locator["index"], a=left, b=right)
    locator.emit(Adjoint(setup))
    location = annotate(locator.finish(), "sparse_location_inplace")

    physical = roe_entry(inputs, **entry_options)
    b = Builder(
        _name("qfvm_sparse_entry", geometry, physical, padding_value),
        {"row": Bits(n), "column": Bits(n), "data": Bits(inputs.fmt.width)},
        resources_for(("geometry", geometry), ("physical", physical)),
        attributes={
            "value_encoding": "signed_fixed_point",
            "value_fraction": inputs.fmt.fraction,
            "matrix": "Hermitian dilation of Roe matrix, identity on padded coordinates",
            "invalid_arithmetic": "entry totalized to zero",
            "correctness": "pending",
        },
    )
    selected = b.local("selected_geometry", Bits(inputs.geometry_width))
    for rank in range(9):
        slot = b.local("slot_" + str(rank), Bits(4))
        geom = b.local("geometry_" + str(rank), Bits(inputs.geometry_width))
        match = b.local("match_" + str(rank), Bits(1))
        for bit in range(4):
            if (rank >> bit) & 1:
                b.x(slot[bit])
        invoke(b, geometry, "geometry", address=fuse(b["column"], slot), data=geom)
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        # 原始 geometry 的 valid 位排除了补齐变量和无效槽位。
        with b.control(fuse(match, geom[inputs.geometry_width - 1])):
            b.xor(geom[: inputs.geometry_width - 1], selected[: inputs.geometry_width - 1])
            b.x(selected[inputs.geometry_width - 1])
        b.call(compare_words(n), a=b["row"], b=geom[:n], flag=match)
        invoke(b, geometry, "geometry", address=fuse(b["column"], slot), data=geom)
    _, _, source, row, col, band, valid = _geometry_refs(selected, inputs)
    value = b.local("computed_value", Bits(inputs.fmt.width))
    status = b.local("arithmetic_status", Bits(2))
    invoke(
        b,
        physical,
        "physical",
        source=source,
        row=row,
        col=col,
        band=band,
        value=value,
        status=status,
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
    entry = annotate(b.finish(), "sparse_entry_xor")
    return SparseAccess(location, entry, n, inputs.fmt.width, 9)
