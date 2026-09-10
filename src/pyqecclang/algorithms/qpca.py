"""QPCA 量子主成分分析：密度矩阵指数化（LMR 协议）+ 相位估计。

实现 Lloyd–Mohseni–Rebentrost 2014（"Quantum principal component analysis"，
Nature Physics 10, 631）的核心原语。密度矩阵指数化利用 SWAP 算符的本征值
结构：对一份 ρ 拷贝与系统做部分交换 ``e^{-iΔt·SWAP}``，丢弃拷贝后系统上
有效通道即 ``e^{-iρΔt}`` 的一阶近似（误差 ``O(Δt²)``）；copies 份拷贝串联
给出总时间 t = copies·Δt。QPCA 对该酉做相位估计，读出 ρ 的谱。

input model：ρ 的拷贝由态制备 oracle 提供（SP/QRAM 三层可绑）；混合态经
density.PurificationAccess.as_state_preparation 适配（环境位不参与交换，
效果等同于取偏迹）。多拷贝意味着对制备 oracle 的多次调用，这是 QPCA 的
资源前提，与 QRAM 假设一并写入属性。
"""

from __future__ import annotations

import math

from pyqecclang.algorithms.contracts import finite_real, positive_integer, require_instance
from pyqecclang.algorithms.fourier import inverse_qft
from pyqecclang.algorithms.operators import _name
from pyqecclang.algorithms.oracles import StatePreparation, invoke, resources_for
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError


def _pauli_pair_rotation(b, first, second, angle, axis):
    """``e^{-i(angle/2)·P⊗P}``，P 由 axis ∈ {x, y, z} 指定；基于奇偶校验的精确分解。"""
    if axis == "x":
        b.h(first)
        b.h(second)
    elif axis == "y":
        # Rx(π/2) = H·Rz(π/2)·H 把 Z 基旋转到 −Y（两项符号相消）。
        for ref in (first, second):
            b.h(ref)
            b.rz(ref, math.pi / 2)
            b.h(ref)
    b.xor(first, second)
    b.rz(second, angle)
    b.xor(first, second)
    if axis == "y":
        for ref in (first, second):
            b.h(ref)
            b.rz(ref, -math.pi / 2)
            b.h(ref)
    elif axis == "x":
        b.h(first)
        b.h(second)


def _partial_swap(b, first, second, angle):
    """``e^{-i·angle·SWAP}``：XX+YY+ZZ = 2·SWAP − I 的精确分解（三轴可交换）。"""
    b.global_phase(-angle / 2)
    for axis in ("z", "x", "y"):
        _pauli_pair_rotation(b, first, second, angle, axis)


def _registers_of(operation):
    return {r.name: r.type for r in operation.module.registers}


def density_matrix_exponentiation(preparation, *, time, copies, swap_width=None, name=None):
    """LMR 密度矩阵指数化：用 copies 份 ρ 拷贝在系统上近似 ``e^{-iρ·time}``。

    Args:
        preparation: ρ 拷贝的态制备（StatePreparation，可为纯化适配结果）。
        time: 总演化时间 t，每步 Δt = t/copies。
        copies: 拷贝数，一阶误差 O(t·Δt) 随之线性下降。
        swap_width: 参与交换的前缀位宽，缺省为制备的整个 target（纯化场景
            应取 ρ 的系统位宽，环境位留在拷贝中不参与交换）。
        name: 覆盖自动生成的模块名。

    Returns:
        Operation: 寄存器 system（swap_width 位）、copies（copies × prep.width 位）。
        系统的输入态由调用方准备；拷贝寄存器从全零由制备 oracle 填充。"""
    require_instance(preparation, StatePreparation, "density_matrix_exponentiation.preparation")
    finite_real(time, "density_matrix_exponentiation.time", minimum=0, strict=True)
    positive_integer(copies, "density_matrix_exponentiation.copies", minimum=1, maximum=63)
    prep_width = preparation.width
    swap_width = prep_width if swap_width is None else swap_width
    positive_integer(
        swap_width, "density_matrix_exponentiation.swap_width", minimum=1, maximum=prep_width
    )
    if preparation.work_width:
        raise ValidationError("拷贝制备需要零 work 的 StatePreparation（纯化适配请先拼接寄存器）")
    step = float(time) / copies
    b = Builder(
        name or _name("dm_exponentiation", preparation.operation, time, copies),
        {"system": Bits(swap_width), "copies": Bits(copies * prep_width)},
        resources_for(("prep", preparation.operation)),
        attributes={
            "algorithm": "density_matrix_exponentiation",
            "copies": copies,
            "step_time": step,
            "error_scaling": "O(time * step_time)",
            "reference": "Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631",
        },
    )
    for k in range(copies):
        copy = b["copies"][k * prep_width : (k + 1) * prep_width]
        invoke(b, preparation.operation, "prep", target=copy, work=copy[:0])
        for bit in range(swap_width):
            _partial_swap(b, b["system"][bit], copy[bit], step)
    return b.finish()


def qpca(preparation, *, precision, step_time, system=None, swap_width=None, name=None):
    """QPCA 主成分分析：对 ``e^{-iρ·step_time}`` 做相位估计，读出 ρ 的谱。

    Args:
        preparation: ρ 拷贝的态制备；共消耗 ``2**precision − 1`` 份拷贝。
        precision: 相位寄存器位数，范围为 1..6（拷贝数随指数增长）。
        step_time: 单步演化时间 Δt；须满足 λ·Δt ≪ 2π 以免读出混叠。
        system: 可选的系统输入态制备（缺省为 ``|0>``）；不同本征态输入
            对应读出不同的本征值峰。
        swap_width: 参与交换的前缀位宽，语义同 density_matrix_exponentiation。
        name: 覆盖自动生成的模块名。

    Returns:
        Operation: 寄存器 system、copies、phase。读出 phase 后用
        eigenvalue_from_phase 解码本征值。"""
    require_instance(preparation, StatePreparation, "qpca.preparation")
    positive_integer(precision, "qpca.precision", minimum=1, maximum=6)
    finite_real(step_time, "qpca.step_time", minimum=0, strict=True)
    prep_width = preparation.width
    swap_width = prep_width if swap_width is None else swap_width
    positive_integer(swap_width, "qpca.swap_width", minimum=1, maximum=prep_width)
    if preparation.work_width:
        raise ValidationError("拷贝制备需要零 work 的 StatePreparation（纯化适配请先拼接寄存器）")
    if system is not None:
        require_instance(system, StatePreparation, "qpca.system")
        if system.width != swap_width:
            raise ValidationError("系统输入态宽度必须等于 swap_width")
    total_copies = (1 << precision) - 1
    resources = [("prep", preparation.operation)]
    if system is not None:
        resources.append(("system", system.operation))
    b = Builder(
        name or _name("qpca", preparation.operation, precision, step_time),
        {
            "system": Bits(swap_width),
            "copies": Bits(total_copies * prep_width),
            "phase": Bits(precision),
        },
        resources_for(*resources),
        attributes={
            "algorithm": "qpca",
            "readout_register": "phase",
            "decoder": "eigenvalue_from_phase",
            "copies": total_copies,
            "step_time": float(step_time),
            "reference": "Lloyd-Mohseni-Rebentrost 2014, Nature Physics 10, 631",
        },
    )
    if system is not None:
        invoke(b, system.operation, "system", target=b["system"], work=b["system"][:0])
    for k in range(total_copies):
        copy = b["copies"][k * prep_width : (k + 1) * prep_width]
        invoke(b, preparation.operation, "prep", target=copy, work=copy[:0])
    b.h(b["phase"])
    cursor = 0
    for bit in range(precision):
        with b.control(b["phase"][bit]):
            for _ in range(1 << bit):
                copy = b["copies"][cursor * prep_width : (cursor + 1) * prep_width]
                for qubit in range(swap_width):
                    _partial_swap(b, b["system"][qubit], copy[qubit], step_time)
                cursor += 1
    invoke(b, inverse_qft(precision), "iqft", target=b["phase"])
    return b.finish()


def eigenvalue_from_phase(value, precision, step_time):
    """把 qpca 的 phase 读出解码为 ρ 的本征值估计。

    酉步 ``e^{-iρΔt}`` 的本征相位为 ``φ = -λΔt/(2π) (mod 1)``；本函数按
    λ ∈ [0, π/Δt) 的分支解码，λ·Δt 超出该范围时发生混叠（调用方责任）。"""
    positive_integer(precision, "eigenvalue_from_phase.precision", maximum=63)
    positive_integer(value, "eigenvalue_from_phase.value", minimum=0, maximum=(1 << precision) - 1)
    finite_real(step_time, "eigenvalue_from_phase.step_time", minimum=0, strict=True)
    if value > 1 << (precision - 1):
        value -= 1 << precision
    return -2.0 * math.pi * value / ((1 << precision) * step_time)
