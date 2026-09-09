"""QFVM 对外提供 CKS 形式的 O_F / O_A，几何表仍只需 O(Ns)。"""

from .builder import Builder
from .ir import Adjoint, Bits, ValidationError, fuse
from .library import _name
from .oracles import SparseAccess, annotate, invoke, resources_for
from .sparse_models import compare_words, value_transposition


def qfvm_sparse_access(inputs, *, padding_value=1.0, **entry_options):
    from .qfvm import _geometry_refs, roe_entry

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
