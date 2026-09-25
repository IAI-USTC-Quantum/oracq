"""Heinrich quantum summation and numerical integration (FO + QRAM input model).

This implements the quantum summation primitive of Heinrich 2002, "Quantum
Summation with an Application to Integration", J. Complexity 18(1), see also
the quantum quadrature rate for function classes in Novak 2001, and assembles
one-dimensional numerical integration from it. The input model is a function
value loader (``|i>|0> -> |i>|f(i)>``), directly compatible with the three
XorDatabase bindings of this repository, abstract, gate and qram.

Readout uses a comparator construction: the threshold register is put into a
uniform superposition and compared with the function value, so the good-state
probability is exactly ``E[v]/2**w``, linear in v with no small-angle
approximation needed; standard amplitude estimation on this marker gives a
quadratic improvement in query complexity over classical Monte Carlo,
``O(1/ε)`` versus ``O(1/ε²)``.
"""

from __future__ import annotations

import math
from collections import namedtuple
from collections.abc import Iterable
from typing import cast

from oracq.algorithms.common.arithmetic import FixedFormat, fixed_arithmetic
from oracq.algorithms.common.estimation import amplitude_from_phase, phase_estimation
from oracq.algorithms.input_model.block_encoding import reflect_zero
from oracq.algorithms.input_model.contracts import (
    finite_real,
    positive_integer,
    require_instance,
)
from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    XorDatabase,
    annotate,
    gate_database,
    invoke,
    qram_database,
    resources_for,
)
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, ValidationError, fuse

LoaderBundle = namedtuple("LoaderBundle", ("database", "memory"))
LoaderBundle.__doc__ = "Summation loader and its simulated memory; memory is empty for the gate binding."


def table_loader(
    values: Iterable[int],
    data_width: int | None = None,
    *,
    backend: str = "gate",
    name: str | None = None,
) -> LoaderBundle:
    """Wrap a nonnegative-integer function table into a summation loader.

    Args:
        values: Nonnegative integer sequence of f(0..N-1), its length zero-padded to a
            power of two.
        data_width: Value word width w, by default the bit length of the maximum value;
            all values must be smaller than 2**w.
        backend: "gate", a truth-table gate implementation, or "qram", a QRAM resource
            whose data table stays out of the IR.
        name: Overrides the automatically generated module name.

    Returns:
        LoaderBundle: database is an XorDatabase with address=index and data=value;
        for backend="qram" memory provides the memory table needed by simulate,
        otherwise it is empty."""
    values = tuple(values)
    if not values:
        raise ValidationError("The function table must not be empty")
    if any(type(v) is not int or v < 0 for v in values):
        raise ValidationError("Function values must be nonnegative integers")
    if data_width is None:
        data_width = max(1, max(values).bit_length())
    positive_integer(data_width, "table_loader.data_width", maximum=64)
    if any(v >= 1 << data_width for v in values):
        raise ValidationError("A function value exceeds the data_width word width")
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
    raise ValidationError("backend must be gate or qram")


def sum_preparation(database: XorDatabase, *, name: str | None = None) -> Operation:
    """State preparation for quantum summation: uniform index plus function value loading
    plus threshold comparison, with the flag=1 probability exactly E[v]/2**w.

    Register layout: target = index(n) | threshold(w) | flag(1), work =
    value(w). The value word stays in work entangled with the index, since the
    comparator readout does not need it cleaned; flag is set by the first of
    the two comparator calls, and after marking on flag the caller should
    invoke this preparation in reverse to restore.

    Args:
        database: Function value loader; the address bit width is the index count and
            the data bit width is the function value bit width.
        name: Name of the generated preparation module; generated automatically from
            the database and widths when omitted.

    Returns:
        Operation: State preparation operation whose flag=1 probability is exactly
        E[v]/2**w.
    """
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


def sum_iterate(database: XorDatabase) -> Operation:
    """Grover iterate for summation: structurally identical to grover_iterate, with the
    marking driven by the flag bit of the preparation target.

    Args:
        database: Function value loader, determining the index count and the function
            value bit width the iterate acts on.

    Returns:
        Operation: Grover iterate operation composed of the preparation, the flag phase
        marking and the zero reflection.
    """
    prep = sum_preparation(database)
    n = cast("int", dict(prep.module.attributes)["index_bits"])
    w = cast("int", dict(prep.module.attributes)["value_bits"])
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


def quantum_sum(
    database: XorDatabase, *, precision: int = 4, name: str | None = None
) -> Operation:
    """Heinrich quantum summation: estimate the mean ``E[f] = (1/N) Σ_i f(i)`` where f
    takes w-bit nonnegative integer values.

    Args:
        database: Function value loader, an XorDatabase with address=index and
            data=value.
        precision: Number of phase register bits, in the range 1..63; the estimation
            error is on the order of O(1/2**precision).

    Returns:
        Operation: Registers target, work and phase. After reading out phase, decode
        the mean with mean_from_phase. The query complexity is O(1/ε), a quadratic
        improvement over the O(1/ε²) of classical Monte Carlo, Heinrich 2002."""
    require_instance(database, XorDatabase, "quantum_sum.database")
    positive_integer(precision, "quantum_sum.precision", maximum=63)
    prep = sum_preparation(database)
    n = cast("int", dict(prep.module.attributes)["index_bits"])
    w = cast("int", dict(prep.module.attributes)["value_bits"])
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


def mean_from_phase(value: int, precision: int, data_width: int) -> float:
    """Decode the phase readout of quantum_sum into the mean estimate ``E[v]`` in integer
    units.

    The good-state probability is p = ``E[v]/2**data_width``, hence
    ``E[v] = amplitude_from_phase * 2**data_width``.

    Args:
        value: Integer count read out from the phase register.
        precision: Number of phase register bits, which must match the quantum_sum
            setting.
        data_width: Function value bit width, in 1..64.

    Returns:
        float: Mean estimate ``E[v]`` in integer value units.
    """
    positive_integer(data_width, "mean_from_phase.data_width", maximum=64)
    return amplitude_from_phase(value, precision) * (1 << data_width)


def heinrich_rate(smoothness: float, dimension: int) -> dict[str, float]:
    """Optimal convergence rate for numerical integration or summation over a function
    class, with error ~ M^{-rate} for M function evaluations.

    Args:
        smoothness: Smoothness parameter s, for example the derivative order of a
            Hölder or Sobolev class; must be positive.
        dimension: Dimension d; must be a positive integer.

    Returns:
        dict: deterministic is s/d, randomized is s/d + 1/2, and quantum is s/d + 1.
        The quantum rate is exactly a quadratic improvement over the randomized
        classical rate, Heinrich 2002 and Novak 2001."""
    finite_real(smoothness, "heinrich_rate.smoothness", minimum=0, strict=True)
    positive_integer(dimension, "heinrich_rate.dimension", minimum=1, maximum=64)
    return {
        "deterministic": smoothness / dimension,
        "randomized": smoothness / dimension + 0.5,
        "quantum": smoothness / dimension + 1.0,
    }


def quantum_integral(
    database: XorDatabase,
    *,
    precision: int = 4,
    interval: float = 1.0,
    name: str | None = None,
) -> Operation:
    """One-dimensional numerical integration: the quantum sum of function values at
    grid points multiplied by the interval length, the composite rectangle rule.

    Function values are quantized to w-bit integers as v/full_scale, with
    full_scale defaulting to 2**data_width − 1. The total error equals the
    discretization error, determined by the grid and the smoothness, see
    heinrich_rate, plus the QAE estimation error. Decoding uses
    integral_from_phase. The interval length must be positive.

    Args:
        database: Loader of function values at grid points, address=index and
            data=value.
        precision: Number of phase register bits, in the range 1..63; the estimation
            error is on the order of O(1/2**precision).
        interval: Integration interval length; must be a positive real number.
        name: Generated module name; generated automatically when omitted.

    Returns:
        Operation: One-dimensional integration summation operation whose decoder is
        integral_from_phase.
    """
    finite_real(interval, "quantum_integral.interval", minimum=0, strict=True)
    operation = quantum_sum(database, precision=precision, name=name)
    attributes = dict(operation.module.attributes)
    attributes.update(
        algorithm="quantum_integral",
        decoder="integral_from_phase",
        interval=float(interval),
    )
    from dataclasses import replace

    from oracq.infrastructure.builder import Operation

    return Operation(
        replace(operation.module, attributes=tuple(sorted(attributes.items()))),
        operation.dependencies,
    )


def integral_from_phase(
    value: int,
    precision: int,
    data_width: int,
    interval: float = 1.0,
    full_scale: float | None = None,
) -> float:
    """Decode the phase readout of quantum_integral into the integral estimate.

    Function values are quantized as v/full_scale, with full_scale defaulting
    to 2**data_width − 1; the integral estimate equals
    ``E[v]/full_scale × interval``.

    Args:
        value: Integer count read out from the phase register.
        precision: Number of phase register bits, which must match the quantum_sum
            setting.
        data_width: Function value bit width, in 1..64.
        interval: Integration interval length; must be a positive real number.
        full_scale: Full-scale value of the function values; defaults to
            2**data_width − 1.

    Returns:
        float: The integral estimate over the interval.
    """
    finite_real(interval, "integral_from_phase.interval", minimum=0, strict=True)
    if full_scale is None:
        full_scale = (1 << data_width) - 1
    finite_real(full_scale, "integral_from_phase.full_scale", minimum=0, strict=True)
    return mean_from_phase(value, precision, data_width) * interval / full_scale
