"""参数化量子电路、MaxCut QAOA 和 VQE 的 Pauli 测量电路。"""

from pyqecclang.algorithms.contracts import finite_real, positive_integer
from pyqecclang.algorithms.interfaces import checked_state_preparation
from pyqecclang.algorithms.operators import _name
from pyqecclang.algorithms.oracles import invoke, resources_for
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


def hardware_efficient_ansatz(width, layers):
    """生成 Ry/Rz 层与相邻 CNOT 组成的参数化电路。

    Args:
        width: target 位宽。
        layers: 参数形状为 [layer][qubit][Ry,Rz]，所有角度以弧度给出。

    Returns:
        Operation: 只有 target。参数优化与重复执行由调用方组织。"""
    positive_integer(width, "ansatz.width", maximum=64)
    layers = tuple(tuple(tuple(pair) for pair in layer) for layer in layers)
    if not layers or any(
        len(layer) != width or any(len(pair) != 2 for pair in layer) for layer in layers
    ):
        raise ValidationError("ansatz 参数需要非空 [layer][qubit][Ry,Rz] 布局")
    for layer in layers:
        for pair in layer:
            for angle in pair:
                finite_real(angle, "ansatz.angle")
    b = Builder(
        _name("hardware_ansatz", width, layers),
        {"target": Bits(width)},
        attributes={"algorithm": "hardware_efficient_ansatz", "layers": len(layers)},
    )
    for layer in layers:
        for bit, (ry, rz) in enumerate(layer):
            b.ry(b["target"][bit], ry)
            b.rz(b["target"][bit], rz)
        for bit in range(width - 1):
            b.xor(b["target"][bit], b["target"][bit + 1])
    return b.finish()


def qaoa_maxcut(width, edges, gammas, betas):
    """生成 MaxCut 的 QAOA cost/mixer 电路。

    Args:
        width: 图的顶点数，也是 target 位宽。
        edges: (u,v,weight) 三元组，权重非负，不接受自环。
        gammas: 各层 cost 演化角。
        betas: 各层 mixer 演化角，长度与 gammas 相同。

    Returns:
        Operation: 从均匀态开始的 QAOA 电路。读取 target 得到一个割的候选位串。

    Cost 为 Σw(1-ZuZv)/2。该函数不运行经典优化器。"""
    positive_integer(width, "qaoa.width", maximum=64)
    edges = tuple(tuple(edge) for edge in edges)
    gammas, betas = tuple(gammas), tuple(betas)
    if not gammas or len(gammas) != len(betas):
        raise ValidationError("QAOA 需要同长且非空的 gamma/beta 列表")
    for angle in gammas + betas:
        finite_real(angle, "qaoa.angle")
    for edge in edges:
        if len(edge) != 3:
            raise ValidationError("MaxCut 每条边为 (u,v,weight)")
        u, v, weight = edge
        positive_integer(u, "qaoa.u", minimum=0, maximum=width - 1)
        positive_integer(v, "qaoa.v", minimum=0, maximum=width - 1)
        finite_real(weight, "qaoa.weight", minimum=0)
        if u == v:
            raise ValidationError("MaxCut 不接受自环")
    b = Builder(
        _name("qaoa_maxcut", width, edges, gammas, betas),
        {"target": Bits(width)},
        attributes={"algorithm": "qaoa_maxcut", "layers": len(gammas)},
    )
    b.h(b["target"])
    for gamma, beta in zip(gammas, betas, strict=True):
        for u, v, weight in edges:
            b.global_phase(-gamma * weight / 2)
            b.xor(b["target"][u], b["target"][v])
            b.rz(b["target"][v], -gamma * weight)
            b.xor(b["target"][u], b["target"][v])
        for bit in range(width):
            b.gate("rx", b["target"][bit], 2 * beta)
    return b.finish()


def pauli_measurement(preparation, word):
    """将制备态旋转到指定 Pauli 测量基。

    Args:
        preparation: 零输入、干净工作区的态制备。
        word: 同宽 I/X/Y/Z 字符串，第一个字符对应最低位。

    Returns:
        Operation: target/work 接口。读取非 I 位的 Z 奇偶性可以估计该 Pauli 字的期望。"""
    prep = checked_state_preparation(preparation)
    if len(word) != prep.width or any(letter not in "IXYZ" for letter in word):
        raise ValidationError("Pauli 测量字必须与态同宽")
    b = Builder(
        _name("pauli_measurement", prep.operation, word),
        {"target": Bits(prep.width), "work": Bits(prep.work_width)},
        resources_for(("prep", prep.operation)),
        attributes={
            "algorithm": "pauli_measurement",
            "pauli_word": word,
            "readout_register": "target",
        },
    )
    invoke(b, prep.operation, "prep", target=b["target"], work=b["work"])
    for bit, letter in enumerate(word):
        if letter == "Y":
            b.gate("phase", b["target"][bit], -1.5707963267948966)
        if letter in "XY":
            b.h(b["target"][bit])
    return b.finish()


def vqe_measurements(preparation, terms):
    """为 VQE 的 Hamiltonian 各项生成测量电路。

    Args:
        preparation: 参数已实例化的 ansatz 或其他态制备。
        terms: (实系数, Pauli 字) 列表。

    Returns:
        tuple: (系数, 测量 Operation) 列表。能量汇总与参数优化在经典侧进行。"""
    circuits = []
    for coefficient, word in terms:
        finite_real(coefficient, "vqe.coefficient")
        circuits.append((coefficient, pauli_measurement(preparation, word)))
    if not circuits:
        raise ValidationError("VQE 需要非空 Hamiltonian 项")
    return tuple(circuits)
