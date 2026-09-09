"""QRAM 原始流场 → Roe 算术 → PTheta → T_L† S T_R → QLSS。"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .arithmetic import FixedFormat
from .builder import Builder
from .ir import Bits, ValidationError
from .library import _name
from .linking import Binding, bind
from .oracles import (
    StatePreparation,
    abstract_database,
    abstract_state_prep,
    annotate,
    invoke,
    qram_database,
    qram_state_prep,
    resources_for,
)
from .roe import ArithmeticBuilder, roe_face


@dataclass(frozen=True)
class RoeQfvmInputs:
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
        return self.cell_width + 3

    @property
    def geometry_width(self):
        return self.width + 4 + self.cell_width + 2 + 2 + 2 + 1


def roe_qfvm_inputs(*, cell_width=2, fmt=None, angle_width=8, prefix="RoeQfvm"):
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
    from .arithmetic import BooleanNetwork

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
    from .qfvm_sparse import qfvm_sparse_access
    from .sparse_models import real_symmetric_sparse_encoding

    if not 0 < padding_value <= amax:
        raise ValidationError("补齐对角值必须为正且不超过元素上界")
    access = qfvm_sparse_access(inputs, padding_value=padding_value, **entry_options)
    return real_symmetric_sparse_encoding(access, inputs.fmt, amax, diagonal_nonnegative=True)


def rhs_qram_preparation(inputs):
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
    from .linking import unresolved

    missing = {r.name for r in unresolved(program)}
    return bind(program, {k: v for k, v in qram_bindings(inputs).items() if k in missing})


def qfvm_memories(inputs, flow, *, amax=8.0):
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
    from .qfvm_sparse import qfvm_sparse_access
    from .qlss import LinearSystem, SparseSystem, SpectralPromise

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
    from .qlss import QLSSProtocol

    if not isinstance(qlss, QLSSProtocol):
        raise ValidationError(
            "请使用声明 input_model 的 QLSSProtocol；裸 callable 不能表达输入模型"
        )
    if spectrum is None:
        raise ValidationError(
            "QFVM 替换 QLSS 需要最小奇异值/范数声明：spectrum=SpectralPromise(...)"
        )
    return qlss.solve(roe_qfvm_problem(inputs, spectrum=spectrum, rhs_norm=rhs_norm, **options))
