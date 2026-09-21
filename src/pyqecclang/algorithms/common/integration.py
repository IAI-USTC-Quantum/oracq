"""Heinrich 量子求和与数值积分（FO + QRAM input model）。

实现 Heinrich 2002（"Quantum Summation with an Application to Integration"，
J. Complexity 18(1)，另见 Novak 2001 的函数类量子求积率）的量子求和原语，
并以此组装一维数值积分。input model 是函数值加载器（``|i>|0> -> |i>|f(i)>``），
与仓库的 XorDatabase 三层绑定（abstract / gate / qram）直接兼容。

读出采用比较器构造：阈值寄存器均匀叠加后与函数值比较，好状态概率恰为
``E[v]/2**w``（对 v 线性，无需小角度近似）；对该标记做标准振幅估计，
查询复杂度相对经典 Monte Carlo 呈二次改进（``O(1/ε)`` 对 ``O(1/ε²)``）。
"""

from __future__ import annotations

import math
from collections import namedtuple

from pyqecclang.algorithms.common.arithmetic import FixedFormat, fixed_arithmetic
from pyqecclang.algorithms.common.estimation import amplitude_from_phase, phase_estimation
from pyqecclang.algorithms.input_model.block_encoding import reflect_zero
from pyqecclang.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from pyqecclang.algorithms.input_model.operators import _name
from pyqecclang.algorithms.input_model.oracles import (
    XorDatabase,
    annotate,
    gate_database,
    invoke,
    qram_database,
    resources_for,
)
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse

LoaderBundle = namedtuple("LoaderBundle", ("database", "memory"))
LoaderBundle.__doc__ = "求和加载器与其模拟内存（gate 绑定时 memory 为空）。"


def table_loader(values, data_width=None, *, backend="gate", name=None):
    """把非负整数函数表包装成求和加载器。

    Args:
        values: f(0..N-1) 的非负整数序列，长度补零到 2 的幂。
        data_width: 值字宽 w，默认取最大值的位长；所有值必须小于 2**w。
        backend: "gate"（真值表门实现）或 "qram"（QRAM 资源，数据表不进 IR）。
        name: 覆盖自动生成的模块名。

    Returns:
        LoaderBundle: database 为 XorDatabase（address=index、data=value）；
        backend="qram" 时 memory 给出 simulate 所需的内存表，否则为空。"""
    values = tuple(values)
    if not values:
        raise ValidationError("函数表不能为空")
    if any(type(v) is not int or v < 0 for v in values):
        raise ValidationError("函数值必须是非负整数")
    if data_width is None:
        data_width = max(1, max(values).bit_length())
    positive_integer(data_width, "table_loader.data_width", maximum=64)
    if any(v >= 1 << data_width for v in values):
        raise ValidationError("函数值越出 data_width 字宽")
    index_bits = max(1, (len(values) - 1).bit_length())
    padded = values + (0,) * ((1 << index_bits) - len(values))
    if backend == "gate":
        return LoaderBundle(
            gate_database(index_bits, data_width, padded, name=name), {}
        )
    if backend == "qram":
        database = qram_database(index_bits, data_width, name=name)
        memory = {
            f"db__{r.name}": dict(enumerate(padded))
            for r in database.operation.module.resources
        }
        return LoaderBundle(database, memory)
    raise ValidationError("backend 必须是 gate 或 qram")


def sum_preparation(database, *, name=None):
    """量子求和的态制备：均匀 index + 函数值加载 + 阈值比较，flag=1 概率恰为 E[v]/2**w。

    寄存器布局：target = index(n) | threshold(w) | flag(1)，work = value(w)。
    value 字与 index 纠缠留在 work（比较器读出不需要复净它）；flag 由比较器
    两次调用中的第一次置位，调用方按 flag 标记后应逆调用本制备以复原。"""
    require_instance(database, XorDatabase, "sum_preparation.database")
    n, w = database.address_width, database.data_width
    compare = fixed_arithmetic("lt", FixedFormat(w, 0, signed=False))
    b = Builder(
        name or _name("sum_prep", database.operation, w),
        {"target": Bits(n + w + 1), "work": Bits(w)},
        resources_for(("db", database.operation), ("lt", compare)),
    )
    index, threshold, flag = b["target"][:n], b["target"][n : n + w], b["target"][n + w :]
    b.h(index)
    b.h(threshold)
    invoke(b, database.operation, "db", address=index, data=b["work"])
    invoke(
        b,
        compare,
        "lt",
        a=threshold,
        b=b["work"],
        out=flag,
        status=b.local("status", Bits(2)),
    )
    return annotate(
        b.finish(),
        "state_prep_isometry",
        zero_input=True,
        clean_work=False,
        implementation="comparator_sum",
        index_bits=n,
        value_bits=w,
    )


def sum_iterate(database):
    """求和的 Grover 迭代：结构同 grover_iterate，标记由制备 target 的 flag 位驱动。"""
    prep = sum_preparation(database)
    n = dict(prep.module.attributes)["index_bits"]
    w = dict(prep.module.attributes)["value_bits"]
    marker_b = Builder(
        _name("sum_mark", n, w), {"target": Bits(n + w + 1)}
    )
    with marker_b.control(marker_b["target"][n + w :], 1):
        marker_b.global_phase(math.pi)
    marker = annotate(marker_b.finish(), "phase_oracle")
    b = Builder(
        _name("sum_iterate", prep, marker),
        {"target": Bits(n + w + 1), "work": Bits(w)},
        resources_for(("prep", prep), ("marker", marker)),
    )
    invoke(b, marker, "marker", target=b["target"])
    with b.adjoint():
        invoke(b, prep, "prep", target=b["target"], work=b["work"])
    reflect_zero(b, fuse(b["target"], b["work"]), positive=True)
    invoke(b, prep, "prep", target=b["target"], work=b["work"])
    return b.finish()


def quantum_sum(database, *, precision=4, name=None):
    """Heinrich 量子求和：估计均值 ``E[f] = (1/N) Σ_i f(i)``，f 取 w 位非负整数值。

    Args:
        database: 函数值加载器（XorDatabase，address=index、data=value）。
        precision: 相位寄存器位数，范围为 1..63；估计误差量级 O(1/2**precision)。

    Returns:
        Operation: 寄存器 target、work、phase。读出 phase 后用
        mean_from_phase 解码均值。查询复杂度 O(1/ε)，相对经典
        Monte Carlo 的 O(1/ε²) 为二次改进（Heinrich 2002）。"""
    require_instance(database, XorDatabase, "quantum_sum.database")
    positive_integer(precision, "quantum_sum.precision", maximum=63)
    prep = sum_preparation(database)
    n = dict(prep.module.attributes)["index_bits"]
    w = dict(prep.module.attributes)["value_bits"]
    iterate = sum_iterate(database)
    qpe = phase_estimation(iterate, precision=precision)
    b = Builder(
        name or _name("quantum_sum", database.operation, precision),
        {"target": Bits(n + w + 1), "work": Bits(w), "phase": Bits(precision)},
        resources_for(("prep", prep), ("qpe", qpe)),
        attributes={
            "algorithm": "quantum_sum",
            "readout_register": "phase",
            "decoder": "mean_from_phase",
            "value_bits": w,
            "index_bits": n,
            "query_complexity": "O(1/epsilon)",
            "classical_query_complexity": "O(1/epsilon**2)",
            "reference": "Heinrich 2002, J. Complexity 18(1)",
        },
    )
    invoke(b, prep, "prep", target=b["target"], work=b["work"])
    invoke(b, qpe, "qpe", target=b["target"], work=b["work"], phase=b["phase"])
    return b.finish()


def mean_from_phase(value, precision, data_width):
    """把 quantum_sum 的 phase 读出解码为均值估计 ``E[v]`` （整数单位）。

    好状态概率 p = ``E[v]/2**data_width``，故 ``E[v] = amplitude_from_phase * 2**data_width``。"""
    positive_integer(data_width, "mean_from_phase.data_width", maximum=64)
    return amplitude_from_phase(value, precision) * (1 << data_width)


def heinrich_rate(smoothness, dimension):
    """函数类数值积分/求和的最优收敛率（误差 ~ M^{-rate}，M 为函数求值次数）。

    Args:
        smoothness: 光滑性参数 s（如 Hölder/Sobolev 类的导数阶），必须为正。
        dimension: 维数 d，必须为正整数。

    Returns:
        dict: deterministic（s/d）、randomized（s/d + 1/2）、quantum（s/d + 1）。
        量子率对随机化经典率恰为二次改进（Heinrich 2002；Novak 2001）。"""
    finite_real(smoothness, "heinrich_rate.smoothness", minimum=0, strict=True)
    positive_integer(dimension, "heinrich_rate.dimension", minimum=1, maximum=64)
    return {
        "deterministic": smoothness / dimension,
        "randomized": smoothness / dimension + 0.5,
        "quantum": smoothness / dimension + 1.0,
    }


def quantum_integral(database, *, precision=4, interval=1.0, name=None):
    """一维数值积分：网格点函数值的量子求和乘以区间长度（复合矩形法则）。

    函数值按 v/full_scale 量化为 w 位整数（默认 full_scale = 2**data_width − 1）。
    总误差 = 离散化误差（由网格/光滑性决定，见 heinrich_rate）+ QAE 估计误差。
    解码用 integral_from_phase。区间长度 interval 必须为正。"""
    finite_real(interval, "quantum_integral.interval", minimum=0, strict=True)
    operation = quantum_sum(database, precision=precision, name=name)
    attributes = dict(operation.module.attributes)
    attributes.update(
        algorithm="quantum_integral",
        decoder="integral_from_phase",
        interval=float(interval),
    )
    from dataclasses import replace

    from pyqecclang.infrastructure.builder import Operation

    return Operation(
        replace(operation.module, attributes=tuple(sorted(attributes.items()))),
        operation.dependencies,
    )


def integral_from_phase(value, precision, data_width, interval=1.0, full_scale=None):
    """把 quantum_integral 的 phase 读出解码为积分估计。

    函数值按 v/full_scale 量化（full_scale 缺省取 2**data_width − 1）；
    积分估计 = ``E[v]/full_scale × interval``。"""
    finite_real(interval, "integral_from_phase.interval", minimum=0, strict=True)
    if full_scale is None:
        full_scale = (1 << data_width) - 1
    finite_real(full_scale, "integral_from_phase.full_scale", minimum=0, strict=True)
    return mean_from_phase(value, precision, data_width) * interval / full_scale
