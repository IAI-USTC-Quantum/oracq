"""Costa QLSS 的 general walk 与相干 filtering 组装原型。

本阶段核对范式、接口和描述可生成性，尚未认证求解精度与成功概率。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ..builder import Builder
from ..combinators import reflect_zero
from ..ir import Bits, ValidationError, fuse
from ..library import BlockEncoding, _name
from ..oracles import StateOracle, StatePreparation, annotate, invoke, resources_for


@dataclass(frozen=True)
class CostaConfig:
    """kappa 为实际编码矩阵 A/alpha 的逆谱界；并非任意尺度下的 cond(A)。"""

    steps: int = 2
    kappa: float = 4.0
    schedule_power: float = 1.5
    filter_degree: int = 2
    filter_attenuation: float = 0.2


@dataclass(frozen=True)
class FilterPlan:
    weights: tuple[float, ...]
    stride: int = 1
    offset: int = 0
    method: str = "explicit"


def dolph_chebyshev_plan(degree=2, attenuation=0.2):
    """以 Laurent 多项式递推构造一个偶数阶 Dolph–Chebyshev LCU。"""
    if type(degree) is not int or degree < 2 or degree % 2:
        raise ValidationError("Dolph–Chebyshev 原型需要正偶数阶")
    if not 0 < attenuation < 1:
        raise ValidationError("attenuation 必须介于 0 和 1 之间")
    beta = math.cosh(math.acosh(1 / attenuation) / degree)
    previous, current = {0: 1.0}, {-1: beta / 2, 1: beta / 2}
    for _ in range(2, degree + 1):
        following = {}
        for power, value in current.items():
            following[power - 1] = following.get(power - 1, 0.0) + beta * value
            following[power + 1] = following.get(power + 1, 0.0) + beta * value
        for power, value in previous.items():
            following[power] = following.get(power, 0.0) - value
        previous, current = current, following
    weights = tuple(
        max(0.0, attenuation * current.get(k, 0.0)) for k in range(-degree, degree + 1, 2)
    )
    return FilterPlan(weights, stride=2, offset=-degree, method="dolph_chebyshev")


def schedule(s, kappa, power=1.5):
    if kappa < 1 or power == 1 or not 0 <= s <= 1:
        raise ValidationError("Costa 调度参数无效")
    if kappa == 1:
        return s
    return kappa / (kappa - 1) * (1 - (1 + s * (kappa ** (power - 1) - 1)) ** (1 / (1 - power)))


def costa_walk(a: BlockEncoding, bprep: StatePreparation, fs: float):
    if a.width != bprep.width or not 0 <= fs <= 1:
        raise ValidationError("Costa 输入形状或调度点不匹配")
    total = a.signal_qubits + bprep.work_width + 4
    b = Builder(
        _name("costa_walk", a.operation, bprep.operation, fs),
        {"target": Bits(a.width), "signal": Bits(total)},
        resources_for(("a", a.operation), ("b", bprep.operation)),
        attributes={"algorithm": "costa_general_walk", "validation_stage": "paradigm", "fs": fs},
    )
    enc = b["signal"][: a.signal_qubits]
    bw = b["signal"][a.signal_qubits : a.signal_qubits + bprep.work_width]
    start = a.signal_qubits + bprep.work_width
    a1, a2, a3, a4 = (b["signal"][start + i] for i in range(4))

    def prep(inverse=False):
        if inverse:
            with b.adjoint():
                invoke(b, bprep.operation, "b", target=b["target"], work=bw)
        else:
            invoke(b, bprep.operation, "b", target=b["target"], work=bw)

    def reflect_rhs_input_zero():
        # U_b 是 target+work 上的酉扩张，投影必须同时要求工作位为零。
        with b.control(b["target"], 0):
            if bw.width:
                with b.control(bw, 0):
                    b.global_phase(math.pi)
            else:
                b.global_phase(math.pi)

    def rotation():
        # R(s) 是反射矩阵，写成 Ry(2 atan2(f,1-f)) Z。
        b.z(a2)
        b.ry(a2, 2 * math.atan2(fs, 1 - fs))

    b.h(a3)
    prep(True)
    b.x(a1)
    with b.control(fuse(a1, a3, a4)):
        reflect_rhs_input_zero()
    b.x(a1)
    prep()
    b.x(a4)
    with b.control(a4):
        rotation()
    b.x(a4)
    with b.control(a4):
        b.h(a2)
    with b.control(fuse(a1, a2)):
        invoke(b, a.operation, "a", target=b["target"], signal=enc)
    with b.control(a2):
        b.x(a1)
    with b.control(a1):
        reflect_zero(b, a2)
    with b.control(fuse(a1, a2)):
        with b.adjoint():
            invoke(b, a.operation, "a", target=b["target"], signal=enc)
    b.x(a4)
    with b.control(a4):
        b.h(a2)
    b.x(a4)
    with b.control(a4):
        rotation()
    b.x(a4)
    prep(True)
    b.x(a1)
    with b.control(fuse(a1, a3, a4)):
        reflect_rhs_input_zero()
    b.x(a1)
    prep()
    b.h(a3)
    reflect_zero(b, fuse(enc, a2, a3), positive=True)
    b.global_phase(math.pi / 2)
    return b.finish()


def unary_weight_preparation(weights):
    weights = tuple(float(v) for v in weights)
    if not weights or any(not math.isfinite(v) or v < 0 for v in weights) or sum(weights) <= 0:
        raise ValidationError("unary PREP 需要非负、非零总和的权重")
    width = len(weights) - 1
    if width < 1:
        raise ValidationError("unary PREP 至少需要两个权重")
    b = Builder(_name("unary_prepare", weights), {"target": Bits(width)})
    for bit in range(width):
        remaining = sum(weights[bit:])
        tail = sum(weights[bit + 1 :])
        angle = 0.0 if remaining == 0 else 2 * math.asin(math.sqrt(tail / remaining))
        if bit:
            with b.control(b["target"][bit - 1]):
                b.ry(b["target"][bit], angle)
        else:
            b.ry(b["target"][bit], angle)
    return annotate(b.finish(), "state_prep_isometry", zero_input=True, encoding="unary_prefix")


def lcu_filter(walk, plan: FilterPlan):
    widths = {r.name: r.type.width for r in walk.module.registers}
    if set(widths) != {"target", "signal"}:
        raise ValidationError("filtering 接收 target/signal walk 接口")
    if type(plan.stride) is not int or plan.stride < 1 or type(plan.offset) is not int:
        raise ValidationError("filter 幂配置无效")
    prep = unary_weight_preparation(plan.weights)
    clock_width = len(plan.weights) - 1
    b = Builder(
        _name("filter", walk, plan.weights, plan.stride, plan.offset),
        {"target": Bits(widths["target"]), "signal": Bits(widths["signal"] + clock_width)},
        resources_for(("walk", walk)),
        attributes={
            "algorithm": "coherent_lcu_filter",
            "filter_method": plan.method,
            "filter_terms": len(plan.weights),
            "stride": plan.stride,
            "offset": plan.offset,
            "validation_stage": "paradigm",
        },
    )
    work, clock = b["signal"][: widths["signal"]], b["signal"][widths["signal"] :]

    def repeat_walk(count):
        with b.repeat(count):
            invoke(b, walk, "walk", target=b["target"], signal=work)

    if plan.offset < 0:
        with b.adjoint():
            repeat_walk(-plan.offset)
    else:
        repeat_walk(plan.offset)
    invoke(b, prep, target=clock)
    for bit in range(clock_width):
        with b.control(clock[bit]):
            repeat_walk(plan.stride)
    with b.adjoint():
        invoke(b, prep, target=clock)
    return b.finish()


def costa_qlss(a: BlockEncoding, bprep: StatePreparation, config=None, *, filtering=None):
    config = CostaConfig() if config is None else config
    if type(config.steps) is not int or config.steps < 1:
        raise ValidationError("Costa steps 必须为正整数")
    walks = tuple(
        costa_walk(a, bprep, schedule((i + 1) / config.steps, config.kappa, config.schedule_power))
        for i in range(config.steps)
    )
    plan = filtering or dolph_chebyshev_plan(config.filter_degree, config.filter_attenuation)
    final_filter = lcu_filter(walks[-1], plan)
    ws = a.signal_qubits + bprep.work_width + 4
    signal_width = ws + len(plan.weights) - 1
    operands = [("prep", bprep.operation), ("filter", final_filter)]
    operands += [(f"walk{i}", walk) for i, walk in enumerate(walks)]
    b = Builder(
        _name("costa_qlss", a.operation, bprep.operation, config, plan),
        {"target": Bits(a.width), "signal": Bits(signal_width)},
        resources_for(*operands),
        attributes={
            "algorithm": "costa_qlss",
            "input_alpha": a.alpha,
            "encoded_inverse_bound": config.kappa,
            "normalization_assumption": "sigma_min(A / alpha) >= 1 / kappa",
            "kernel_status": "prototype; initial walk eigenstate and readout channel unverified",
            "validation_stage": "paradigm",
            "steps": config.steps,
            "filtering": plan.method,
            "success_condition": "signal == 0; probability and solution accuracy unverified",
        },
    )
    bw = b["signal"][a.signal_qubits : a.signal_qubits + bprep.work_width]
    invoke(b, bprep.operation, "prep", target=b["target"], work=bw)
    for i, walk in enumerate(walks):
        invoke(b, walk, f"walk{i}", target=b["target"], signal=b["signal"][:ws])
    invoke(b, final_filter, "filter", target=b["target"], signal=b["signal"])
    return StateOracle(b.finish())


def make_costa_qlss(config=None):
    """声明 BE 输入；问题层自动按 alpha/sigma_min 推导实际调度参数。"""
    from dataclasses import replace

    from ..qlss import QLSSProtocol

    config = CostaConfig() if config is None else config

    def kernel(system):
        effective = replace(config, kappa=system.inverse_norm_bound)
        return costa_qlss(system.encoding, system.rhs, effective)

    return QLSSProtocol(
        "costa_general_walk", "block_encoding", kernel, legacy=lambda a, b: costa_qlss(a, b, config)
    )
