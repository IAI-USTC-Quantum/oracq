"""Quantum building blocks of the QCNN (arXiv:1911.01117 §5.1).

``qcnn_vector_prep`` implements the QRAM row/column loading of
Eq. (13)-(15): the index register p selects a row and the target register
prepares the normalized row vector state ``|A_p⟩``. Rotation angles and signs
are stored as runtime QRAM banks — an angle tree (RY half-angle convention,
addressing (1<<depth)-1+prefix) plus sign kickback, a mechanism consistent
with the QFVM sign residual tree; swapping data only swaps the memory tables
and leaves the circuit unchanged.

``qcnn_inner_product`` implements Eq. (16)-(19): p and q are put into uniform
superposition and a Hadamard test on flag encodes the inner product into an
amplitude. The probability of measuring (p, q, flag) satisfies
P0(p,q) = (1+<A_p|F_q>)/(2 H'W'D') (Eq. 20).

The semantic validation of both building blocks against a classical mirror
with quantized angles lives in tests/core/test_qcnn.py.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from oracq.algorithms.input_model.operators import _name
from oracq.algorithms.input_model.oracles import (
    XorDatabase,
    invoke,
    qram_database,
    resources_for,
)
from oracq.algorithms.qml.qcnn import ConvSpec
from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Ref, ValidationError, fuse


def _log2_exact(value: int) -> int:
    """Validate that ``value`` is a power of two and return its base-2 logarithm."""
    if value < 1 or value & (value - 1):
        raise ValidationError("QCNN building blocks require the row and column count to be a power of two")
    return value.bit_length() - 1


def vector_angle_tables(
    rows: Sequence[Sequence[float]], angle_width: int
) -> tuple[dict[int, int], dict[int, int]]:
    """Compile a group of equal-length row vectors into an angle-tree bank and a sign bank (runtime memory tables).

    Angles are derived from the squared partial sums under the RY half-angle
    convention: the RY parameter is 2*atan2(sqrt(S_right), sqrt(S_left)); tree
    nodes are addressed (1<<depth)-1+prefix (matching the addressing of
    qram_state_prep). Returns (angles, signs): angles is keyed by
    fuse(index, node) and signs by fuse(index, leaf).

    Args:
        rows: Sequence of equal-length nonzero-norm row vectors; each vector
            length must be a power of two.
        angle_width: Quantization bit width of an angle word; angles are
            quantized in units of 2π/2^angle_width radians.

    Returns:
        tuple[dict[int, int], dict[int, int]]: Runtime memory tables of the
        angle-tree bank and the sign bank.
    """
    length = len(rows[0])
    if length < 1 or length & (length - 1):
        raise ValidationError("the vector length must be a power of two")
    width = length.bit_length() - 1
    angles: dict[int, int] = {}
    signs: dict[int, int] = {}
    scale = (1 << angle_width) / (2 * math.pi)
    for index, row in enumerate(rows):
        if any(len(r) != length for r in rows):
            raise ValidationError("all rows must have the same length")
        squared = [value * value for value in row]
        if sum(squared) == 0:
            raise ValidationError("the row vector has zero norm")
        for depth in range(width):
            bit = width - depth - 1
            for prefix in range(1 << depth):
                start = prefix << (bit + 1)
                left = sum(squared[start : start + (1 << bit)])
                right = sum(squared[start + (1 << bit) : start + (1 << (bit + 1))])
                node = (1 << depth) - 1 + prefix
                angle = 2 * math.atan2(math.sqrt(right), math.sqrt(left))
                angles[(index << width) | node] = round(angle * scale) % (1 << angle_width)
        for leaf, value in enumerate(row):
            signs[(index << width) | leaf] = int(value < 0)
    return angles, signs


def _constant_register(
    builder: Builder, name: str, width: int, value: int
) -> Ref:
    """A constant register prepared with X gates (starting from |0>; the caller handles the symmetric restoration)."""
    register = builder.local(name, Bits(width))
    for bit in range(width):
        if (value >> bit) & 1:
            builder.x(register[bit])
    return register


def _prep_node_rotation(
    builder: Builder,
    angles: XorDatabase,
    node_register: Ref,
    bit: int,
    angle_width: int,
) -> None:
    """The three steps on a single tree node: angle lookup, controlled RY synthesis, and uncomputing re-lookup.

    The bank key is (index<<width)|node: index occupies the high bits, so the
    address is fuse(node_register, index).
    """
    invoke(
        builder,
        angles.operation,
        "angles",
        address=fuse(node_register, builder["index"]),
        data=builder["work"],
    )
    for k in range(angle_width):
        with builder.control(builder["work"][k], 1):
            builder.ry(
                builder["target"][bit], 2 * math.pi * (1 << k) / (1 << angle_width)
            )
    invoke(
        builder,
        angles.operation,
        "angles",
        address=fuse(node_register, builder["index"]),
        data=builder["work"],
    )


def qcnn_vector_prep(
    count: int, vector_width: int, angle_width: int, *, name: str | None = None
) -> Operation:
    """Eq. (15): write the normalized vector state into the target register according to the index (QRAM data driven).

    Registers: index (log2 count bits), target (vector_width bits), and work
    (angle_width bits). Resources: angles (address fuse(index, node), data
    width angle_width) and signs (address fuse(index, leaf), data width 1).
    Angle words are synthesized bit by bit into controlled RY rotations:
    2*pi*2^k/2^angle_width. Signs are written via Z kickback and only a branch
    phase remains after the query is uncomputed. Behavior for index >= count
    is undefined; the caller must ensure the superposition only covers valid
    rows.

    Args:
        count: Number of vectors, must be a power of two; determines the index
            register width.
        vector_width: Quantum bit width of a single vector (the vector length
            is a power of two).
        angle_width: Word width of the angle bank; rotations are synthesized
            from it by bit-weighted summation.
        name: Name of the generated module; generated from the parameters by
            default.

    Returns:
        Operation: A preparation operation that selects a row by index and
        writes the normalized vector state into target.
    """
    index_width = _log2_exact(count)
    angles = qram_database(index_width + vector_width, angle_width)
    signs = qram_database(index_width + vector_width, 1)
    b = Builder(
        name or _name("qcnn_vector_prep", count, vector_width, angle_width),
        {
            "index": Bits(index_width),
            "target": Bits(vector_width),
            "work": Bits(angle_width),
        },
        resources_for(("angles", angles.operation), ("signs", signs.operation)),
        attributes={
            "algorithm": "qcnn_qram_row_preparation",
            "implementation": "angle_tree_with_sign_kickback",
            "qram_queries": 2 * vector_width,
            "correctness": "pending",
        },
    )
    for depth in range(vector_width):
        bit = vector_width - depth - 1
        for prefix in range(1 << depth):
            node = (1 << depth) - 1 + prefix
            node_register = _constant_register(b, f"node_{depth}_{prefix}", vector_width, node)
            if depth:
                with b.control(b["target"][bit + 1 :], prefix):
                    _prep_node_rotation(b, angles, node_register, bit, angle_width)
            else:
                _prep_node_rotation(b, angles, node_register, bit, angle_width)
            for bit_ in range(vector_width):
                if (node >> bit_) & 1:
                    b.x(node_register[bit_])
    invoke(
        b,
        signs.operation,
        "signs",
        address=fuse(b["target"], b["index"]),
        data=b["work"][:1],
    )
    b.z(b["work"][0])
    invoke(
        b,
        signs.operation,
        "signs",
        address=fuse(b["target"], b["index"]),
        data=b["work"][:1],
    )
    return b.finish()


def qcnn_inner_product(
    row_count: int,
    column_count: int,
    vector_width: int,
    angle_width: int,
    *,
    name: str | None = None,
) -> Operation:
    """Eq. (16)-(19): Hadamard-type inner product circuit.

    Registers: p and q (row/column indices), flag (1 bit), vec (vector_width
    bits), and work (angle_width + 1 bits). After p and q are put into uniform
    superposition, the flag=0 branch loads a row vector and the flag=1 branch
    loads a column vector; the final Hadamard encodes (1+<A_p|F_q>)/2 into the
    flag=0 amplitude. Resources: one angle and sign bank set each for rows and
    columns (resource name prefixes row and col).

    Args:
        row_count: Number of row vectors, must be a power of two.
        column_count: Number of column vectors, must be a power of two.
        vector_width: Quantum bit width of the row/column vectors.
        angle_width: Word width of the angle bank, matching the preparation
            operation.
        name: Name of the generated module; generated from the parameters by
            default.

    Returns:
        Operation: An inner product estimation operation that encodes
        (1+<A_p|F_q>)/2 into the flag amplitude.
    """
    row_prep = qcnn_vector_prep(row_count, vector_width, angle_width, name="qcnn_row_prep")
    column_prep = qcnn_vector_prep(
        column_count, vector_width, angle_width, name="qcnn_col_prep"
    )
    b = Builder(
        name or _name("qcnn_inner_product", row_count, column_count, vector_width, angle_width),
        {
            "p": Bits(_log2_exact(row_count)),
            "q": Bits(_log2_exact(column_count)),
            "flag": Bits(1),
            "vec": Bits(vector_width),
            "work": Bits(angle_width),
        },
        resources_for(
            ("row", row_prep),
            ("col", column_prep),
        ),
        attributes={
            "algorithm": "qcnn_hadamard_inner_product",
            "equation": "P(p,q,flag=0) = (1+<A_p|F_q>)/(2*rows*cols)",
            "correctness": "pending",
        },
    )
    for register in (b["p"], b["q"], b["flag"]):
        b.h(register)
    with b.control(b["flag"], 0):
        invoke(
            b,
            row_prep,
            "row",
            index=b["p"],
            target=b["vec"],
            work=b["work"],
        )
    with b.control(b["flag"], 1):
        invoke(
            b,
            column_prep,
            "col",
            index=b["q"],
            target=b["vec"],
            work=b["work"],
        )
    b.h(b["flag"])
    return b.finish()

def quantized_prepared_state(vector: Sequence[float], angle_width: int) -> list[float]:
    """Angle-quantization mirror: the (quantized) normalized vector prepared from a row/column.

    Consistent with the circuit semantics: layer-by-layer rotation, each
    layer's angle synthesized bit by bit from the quantized angle word in the
    bank; signs applied per the sign bank. Returns a list of reals.

    Args:
        vector: Real vector to prepare, of power-of-two length and nonzero
            norm.
        angle_width: Quantization bit width of an angle word, matching the
            circuit bank.

    Returns:
        list[float]: The quantized normalized vector matching the circuit
        output.
    """
    length = len(vector)
    if length & (length - 1):
        raise ValidationError("the vector length must be a power of two")
    width = length.bit_length() - 1
    step = 2 * math.pi / (1 << angle_width)
    squared = [value * value for value in vector]
    if sum(squared) == 0:
        raise ValidationError("the row vector has zero norm")
    states = {0: 1.0}
    for depth in range(width):
        bit = width - depth - 1
        nxt: dict[int, float] = {}
        for prefix, amplitude in states.items():
            start = prefix << (bit + 1)
            left = sum(squared[start : start + (1 << bit)])
            right = sum(squared[start + (1 << bit) : start + (1 << (bit + 1))])
            angle = 2 * math.atan2(math.sqrt(right), math.sqrt(left))
            word = round(angle / step) % (1 << angle_width)
            theta = sum(
                2 * math.pi * (1 << k) / (1 << angle_width)
                for k in range(angle_width)
                if (word >> k) & 1
            )
            low, high = prefix * 2, prefix * 2 + 1
            nxt[low] = nxt.get(low, 0.0) + amplitude * math.cos(theta / 2)
            nxt[high] = nxt.get(high, 0.0) + amplitude * math.sin(theta / 2)
        states = nxt
    magnitude = [states[index] for index in range(length)]
    norm = math.sqrt(sum(value * value for value in magnitude))
    return [m / norm * (1.0 if vector[index] >= 0 else -1.0) for index, m in enumerate(magnitude)]


def qcnn_sampled_layer(
    x: Sequence[float],
    kernel: Sequence[float],
    spec: ConvSpec,
    *,
    samples: int,
    eta: float,
    seed: int = 0,
) -> tuple[list[list[float]], dict[str, int]]:
    """Sampling driver for Algorithm 1 of the paper (small scale, end to end).

    The quantum semantics split into two layers: row/column preparation and
    the Hadamard inner product are already validated by circuits
    (qcnn_inner_product, Eq. 16-20); conditional rotation, amplitude
    amplification, and l_inf tomography (Eq. 28-38) are modeled here by their
    exact distributions — the sampling probability is f(Y_pq)^2/sum f^2
    (Eq. 34), each sample yields the triple (p, q, f(Y_pq)), and pixels below
    eta are treated as unsampled and set to zero (Eq. 36). Sampled values are
    mapped into pooling regions per Eq. 39 and turned into the output tensor
    through the QRAM online overwrite rule (max keeps the high value / average
    spreads it). Returns (pooled_rows, stats).

    Args:
        x: Input tensor data stored flat in row-major order.
        kernel: Flat sequence of convolution kernel coefficients.
        spec: Convolution layer specification, determining the unrolling, the
            activation cap, and the pooling kind.
        samples: Number of samples, a positive integer.
        eta: Pixel retention threshold; sampled values below it are treated as
            unsampled and set to zero.
        seed: Random number generator seed.

    Returns:
        tuple[list[list[float]], dict[str, int]]: The pooled output matrix and
        sampling statistics (sampled/kept counts).
    """
    import random

    from oracq.algorithms.qml.qcnn import (
        QCNNQRAM,
        cap_relu,
        im2col,
        kernel_columns,
    )

    rows = im2col(x, spec)
    columns = kernel_columns(kernel, spec)
    if len(rows) & (len(rows) - 1) or len(columns) & (len(columns) - 1):
        raise ValidationError("the sampling driver requires the row and column counts to be powers of two")
    qram = QCNNQRAM(rows)
    prepared_rows = [quantized_prepared_state(row, 10) for row in rows]
    prepared_columns = [quantized_prepared_state(column, 10) for column in columns]
    values: dict[tuple[int, int], float] = {}
    for p, a in enumerate(prepared_rows):
        for q, f in enumerate(prepared_columns):
            overlap = sum(u * v for u, v in zip(a, f, strict=True))
            inner = overlap * qram.norm(p) * _column_norm(columns[q])
            values[(p, q)] = cap_relu(inner, spec.cap)
    weight = {key: value * value for key, value in values.items()}
    total = sum(weight.values())
    if total == 0:
        return _zero_pooled(spec), {"sampled": 0, "kept": 0}
    rng = random.Random(seed)
    oh, ow, _ = spec.output_shape
    ph, pw, d = spec.pooled_shape
    pooled = [[0.0] * d for _ in range(ph * pw)]
    counts = [[0] * d for _ in range(ph * pw)]
    kept = 0
    for _ in range(samples):
        pick, acc = rng.random() * total, 0.0
        for key, w in weight.items():
            acc += w
            if acc >= pick:
                p, q = key
                break
        value = values[(p, q)]
        if value < eta:
            continue
        kept += 1
        i, j = p // ow, p % ow
        index = (i // spec.pool) * pw + (j // spec.pool)
        if spec.pool_kind == "max":
            pooled[index][q] = max(pooled[index][q], value)
        else:
            pooled[index][q] = (pooled[index][q] * counts[index][q] + value) / (counts[index][q] + 1)
        counts[index][q] += 1
    stats = {"sampled": samples, "kept": kept}
    return pooled, stats


def _column_norm(column: Sequence[float]) -> float:
    """Return the Euclidean norm of a column vector."""
    return math.sqrt(sum(value * value for value in column))


def _zero_pooled(spec: ConvSpec) -> list[list[float]]:
    """An all-zero pooled output tensor matching ``spec.pooled_shape``."""
    ph, pw, d = spec.pooled_shape
    return [[0.0] * d for _ in range(ph * pw)]
