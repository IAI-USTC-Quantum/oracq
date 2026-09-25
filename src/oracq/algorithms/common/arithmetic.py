"Reversible fixed-point arithmetic: Boolean SSA to compute/XOR/uncompute, plus native implementations of the same graph."

from __future__ import annotations

import hashlib
import json
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from functools import lru_cache
from typing import TYPE_CHECKING, cast

from oracq.infrastructure.builder import Builder, Operation
from oracq.infrastructure.ir import Bits, Program, Ref, ValidationError, fuse

if TYPE_CHECKING:
    from types import ModuleType

    from oracq.infrastructure.native import NativeContext, NativeRegistry


@dataclass(frozen=True)
class FixedFormat:
    """Describes the fixed-point number format used by reversible arithmetic.

    The low ``fraction`` bits of the ``width``-bit word hold the fraction;
    when ``signed`` is True the most significant bit is the two's-complement
    sign bit. Values and bit patterns are converted through ``encode`` and
    ``decode``.

    Attributes:
        width: Total bit width, including the sign bit if present; the valid range is 2..64.
        fraction: Number of fraction bits, which must not occupy the sign bit.
        signed: Whether the most significant bit is a two's-complement sign bit.

    Raises:
        ValidationError: ``width`` falls outside 2..64, or ``fraction`` occupies the sign bit.
    """
    width: int = 8
    fraction: int = 3
    signed: bool = True

    def __post_init__(self) -> None:
        """Validate the static constraints on the width and fraction combination."""
        if not 2 <= self.width <= 64 or not 0 <= self.fraction < self.width - int(self.signed):
            raise ValidationError("Fixed-point format requires 2≤width≤64 and fraction bits that"
                                  " do not occupy the sign bit")

    def encode(self, value: float) -> int:
        """Encode a fixed-point value as a ``width``-bit integer.

        The value is first multiplied by ``2**fraction`` and truncated toward
        zero, then masked to ``width`` bits.

        Args:
            value: The fixed-point value to encode.

        Returns:
            int: The integer bit pattern after wrapping modulo ``2**width``.
        """
        return int(value * (1 << self.fraction)) & ((1 << self.width) - 1)

    def decode(self, value: int) -> float:
        """Decode a ``width``-bit integer bit pattern back to a fixed-point value.

        Args:
            value: The integer bit pattern; higher bits beyond the word width are masked off first.

        Returns:
            float: The value after interpreting the sign bit as two's complement in the signed
            format, then dividing by ``2**fraction``.
        """
        value &= (1 << self.width) - 1
        if self.signed and value >> (self.width - 1):
            value -= 1 << self.width
        return value / (1 << self.fraction)


class BooleanNetwork:
    """0/1 are constants; every other index is an input or an immutable Boolean value."""

    def __init__(self) -> None:
        """Create an empty network containing only the two constants 0 and 1."""
        self.nodes: list[tuple[str | int, ...]] = [("const", 0), ("const", 1)]
        self.cache: dict[tuple[str | int, ...], int] = {}
        self.inputs: dict[str, list[int]] = {}
        self.outputs: dict[str, list[int]] = {}

    def node(self, key: Iterable[str | int]) -> int:
        """Return the node index for a key, creating a new node when absent.

        Hash-consing on tuple keys: identical keys always return the same index.

        Args:
            key: Iterable key describing a constant, input, or operation.

        Returns:
            int: The node index of the key.
        """
        key = tuple(key)
        if key not in self.cache:
            self.cache[key] = len(self.nodes)
            self.nodes.append(key)
        return self.cache[key]

    def input(self, name: str, width: int) -> list[int]:
        """Declare an input port and register its per-bit nodes.

        Args:
            name: Port name.
            width: Port bit width.

        Returns:
            list: List of node indices, indexed by bit position, least significant bit first.
        """
        bits = [self.node(("input", name, i)) for i in range(width)]
        self.inputs[name] = bits
        return bits

    def inv(self, a: int) -> int:
        """Invert a node, with constant folding and double-negation simplification.

        Args:
            a: Node index.

        Returns:
            int: Node index of the inverted result.
        """
        if a < 2:
            return 1 - a
        if self.nodes[a][0] == "not":
            return cast("int", self.nodes[a][1])
        return self.node(("not", a))

    def xor(self, a: int, b: int) -> int:
        """XOR two nodes, with the constant identity and ``a xor a = 0`` simplifications.

        Args:
            a: Node index.
            b: Node index.

        Returns:
            int: Node index of the XOR result.
        """
        if a == b:
            return 0
        if a == 0 or b == 0:
            return a or b
        if a == 1 or b == 1:
            return self.inv(b if a == 1 else a)
        return self.node(("xor", *sorted((a, b))))

    def and_(self, a: int, b: int) -> int:
        """Bitwise AND of two nodes, with constant absorption and idempotence simplifications.

        Args:
            a: Node index.
            b: Node index.

        Returns:
            int: Node index of the AND result.
        """
        if a == b:
            return a
        if a == 0 or b == 0:
            return 0
        if a == 1 or b == 1:
            return b if a == 1 else a
        return self.node(("and", *sorted((a, b))))

    def or_(self, a: int, b: int) -> int:
        """Bitwise OR of two nodes, composed through ``xor`` and ``and_``.

        Args:
            a: Node index.
            b: Node index.

        Returns:
            int: Node index of the OR result.
        """
        return self.xor(self.xor(a, b), self.and_(a, b))

    def any(self, bits: Iterable[int]) -> int:
        """Reduce a group of nodes to their bitwise OR.

        Args:
            bits: Iterable of node indices.

        Returns:
            int: Node index of the reduced result; an empty input yields constant 0.
        """
        result = 0
        for bit in bits:
            result = self.or_(result, bit)
        return result

    def mux(self, select: int, yes: list[int], no: list[int]) -> list[int]:
        """Select bit by bit between two equal-width signals according to a select bit.

        Args:
            select: Node index of the select bit.
            yes: Bit list to output when the select bit is 1.
            no: Bit list to output when the select bit is 0; must match the width of ``yes``.

        Returns:
            list: Node indices of the bit-by-bit selection result.
        """
        return [
            self.xor(b, self.and_(select, self.xor(a, b))) for a, b in zip(yes, no, strict=True)
        ]

    @staticmethod
    def const(value: int, width: int) -> list[int]:
        """Expand a nonnegative integer into a ``width``-bit constant bit list.

        Args:
            value: The nonnegative integer constant.
            width: Bit width.

        Returns:
            list: Each bit is 0 or 1, least significant bit first.
        """
        return [(value >> i) & 1 for i in range(width)]

    @staticmethod
    def resize(bits: list[int], width: int) -> list[int]:
        """Truncate a bit list or zero-pad its high end to the given width.

        Args:
            bits: The original bit list.
            width: The target bit width.

        Returns:
            list: A new list whose length is exactly ``width``.
        """
        return (list(bits) + [0] * width)[:width]

    def add(self, a: list[int], b: list[int], carry: int = 0) -> tuple[list[int], int]:
        """Ripple-carry addition.

        Args:
            a: Bit list of one addend.
            b: Bit list of the other addend, the same width as ``a``.
            carry: Node index of the initial carry entering the least significant bit.

        Returns:
            tuple: ``(result, carry)``, the bit-by-bit sum and the carry out of the most
            significant bit.
        """
        result = []
        for x, y in zip(a, b, strict=True):
            p = self.xor(x, y)
            result.append(self.xor(p, carry))
            carry = self.xor(self.and_(x, y), self.and_(p, carry))
        return result, carry

    def neg(self, a: list[int]) -> list[int]:
        """Two's-complement negation: invert every bit, then add one at the least significant bit.

        Args:
            a: Bit list of the operand.

        Returns:
            list: Bit list of ``-a``.
        """
        return self.add([self.inv(x) for x in a], [0] * len(a), 1)[0]

    def sub(self, a: list[int], b: list[int]) -> list[int]:
        """Subtraction, implemented as ``a + ~b + 1``.

        Args:
            a: Bit list of the minuend.
            b: Bit list of the subtrahend, the same width as ``a``.

        Returns:
            list: Bit list of ``a - b``.
        """
        return self.add(a, [self.inv(x) for x in b], 1)[0]

    def lt(self, a: list[int], b: list[int]) -> int:
        """Unsigned comparison ``a < b``.

        Args:
            a: Bit list of the left operand.
            b: Bit list of the right operand, the same width as ``a``.

        Returns:
            int: Node index of the comparison result, truth value 1.
        """
        return self.inv(self.add(a, [self.inv(x) for x in b], 1)[1])

    def abs(self, a: list[int], signed: bool) -> list[int]:
        """Return the absolute value of the operand.

        Args:
            a: Bit list of the operand.
            signed: When True, select between ``a`` and ``-a`` by the sign bit; when False,
                return ``a`` unchanged.

        Returns:
            list: Bit list of the absolute value.
        """
        return self.mux(a[-1], self.neg(a), a) if signed else a

    def mul(self, a: list[int], b: list[int], width: int) -> list[int]:
        """Shift-and-accumulate multiplication, with the product truncated to ``width`` bits.

        Args:
            a: Bit list of the multiplicand.
            b: Bit list of the multiplier.
            width: Result bit width.

        Returns:
            list: Bit list of the low ``width`` bits of the product.
        """
        result = [0] * width
        for i, bit in enumerate(b[:width]):
            row = [0] * i + [self.and_(bit, x) for x in a[: width - i]]
            result = self.add(result, self.resize(row, width))[0]
        return result

    def div(self, a: list[int], b: list[int]) -> list[int]:
        """Restoring division, returning only the quotient.

        The dividend bits are shifted into the partial remainder from the most
        significant end; whenever the remainder is at least the divisor, the
        divisor is subtracted and the quotient bit is set.

        Args:
            a: Bit list of the dividend.
            b: Bit list of the divisor.

        Returns:
            list: Bit list of the quotient, with the same width as ``a``; the remainder is
            discarded.
        """
        size = max(len(a), len(b)) + 1
        denominator = self.resize(b, size)
        remainder = [0] * size
        result = [0] * len(a)
        for i in reversed(range(len(a))):
            remainder = [a[i]] + remainder[:-1]
            accept = self.inv(self.lt(remainder, denominator))
            remainder = self.mux(accept, self.sub(remainder, denominator), remainder)
            result[i] = accept
        return result

    def sqrt(self, a: list[int]) -> list[int]:
        """Digit-by-digit restoring square root.

        When the input has an odd number of bits, it is first zero-padded at
        the high end to an even count; the radicand bits then enter in
        two-bit groups from the high end.

        Args:
            a: Bit list of the radicand.

        Returns:
            list: Bit list of the square root, with half the padded input bit count.
        """
        a = self.resize(a, len(a) + len(a) % 2)
        n = len(a) // 2
        size = n + 2
        remainder, root = [0] * size, [0] * size
        for i in reversed(range(n)):
            remainder = [a[2 * i], a[2 * i + 1]] + remainder[:-2]
            trial = [1, 0] + root[:-2]
            accept = self.inv(self.lt(remainder, trial))
            remainder = self.mux(accept, self.sub(remainder, trial), remainder)
            root = [accept] + root[:-1]
        return root[:n]

    def evaluate(self, **inputs: int) -> dict[str, int]:
        """Evaluate the whole network on the classical side.

        Args:
            **inputs: Mapping from port names to nonnegative integers, read bit by bit
                into the input ports.

        Returns:
            dict: Each output port name mapped to its integer evaluation result.
        """
        values = [0, 1]
        for op, *args in self.nodes[2:]:
            if op == "input":
                name, bit = args
                v = (inputs[cast("str", name)] >> cast("int", bit)) & 1
            elif op == "not":
                v = 1 ^ values[cast("int", args[0])]
            elif op == "xor":
                v = values[cast("int", args[0])] ^ values[cast("int", args[1])]
            else:
                v = values[cast("int", args[0])] & values[cast("int", args[1])]
            values.append(v)
        return {
            name: sum(values[v] << i for i, v in enumerate(bits))
            for name, bits in self.outputs.items()
        }

    def payload(self) -> str:
        """Serialize the network into a compact JSON string.

        Returns:
            str: The payload containing nodes, inputs and output ports, usable as a cache
            key or module attribute.
        """
        return json.dumps(
            dict(nodes=self.nodes, inputs=self.inputs, outputs=self.outputs), separators=(",", ":")
        )

    @classmethod
    def from_payload(cls, value: str) -> BooleanNetwork:
        """Rebuild a network from a JSON payload, validating the payload fully as
        untrusted input.

        The validation covers the constant prefix, port naming and mutual
        exclusion, bit widths in 1..64, bit reference ranges, the sequential
        directed-acyclic property, and the consistency of the input port
        mapping.

        Args:
            value: The JSON string produced by ``payload``.

        Returns:
            BooleanNetwork: The network rebuilt after the validation passes.

        Raises:
            ValidationError: Any structural or consistency check of the payload fails.
        """
        data = json.loads(value)
        net = cls()
        net.nodes = [tuple(x) for x in data["nodes"]]
        net.inputs, net.outputs = data["inputs"], data["outputs"]
        from oracq.infrastructure.validation import name

        if net.nodes[:2] != [("const", 0), ("const", 1)] or set(net.inputs) & set(net.outputs):
            raise ValidationError("Invalid Boolean network constants or ports")
        for port, bits in {**net.inputs, **net.outputs}.items():
            name(port)
            if not isinstance(bits, list) or not 1 <= len(bits) <= 64:
                raise ValidationError("Invalid Boolean network port bit width")
            if any(type(v) is not int or not 0 <= v < len(net.nodes) for v in bits):
                raise ValidationError("Boolean network bit reference out of range")
        for index, node in enumerate(net.nodes[2:], 2):
            op, *args = node
            if op == "input":
                if len(args) != 2 or args[0] not in net.inputs or type(args[1]) is not int:
                    raise ValidationError("Invalid Boolean network input node")
                if (
                    not 0 <= args[1] < len(net.inputs[args[0]])
                    or net.inputs[args[0]][args[1]] != index
                ):
                    raise ValidationError("Inconsistent Boolean network input mapping")
            elif op not in {"not", "and", "xor"} or len(args) != (1 if op == "not" else 2):
                raise ValidationError("Unknown Boolean network node")
            elif any(type(v) is not int or not 0 <= v < index for v in args):
                raise ValidationError("The Boolean network is not a sequential directed acyclic"
                                      " graph")
        for port, bits in net.inputs.items():
            if any(net.nodes[v] != ("input", port, i) for i, v in enumerate(bits)):
                raise ValidationError("Incomplete Boolean network input port")
        return net

    def operation(
        self,
        name: str | None = None,
        *,
        attributes: Mapping[str, str | int | float | bool] | None = None,
    ) -> Operation:
        """Compile the network into a reversible quantum operation in compute/copy/uncompute
        form.

        Every operation node is first computed forward into a private
        ``ssa_``-prefixed bank grouped in 64-bit words (``not`` and ``xor``
        via X/XOR copies, ``and`` via a controlled X); after the output bits
        are XOR-copied into the public registers, the adjoint of the forward
        compute frame is emitted to clean all banks, satisfying the
        ``zero_in_zero_out`` workspace contract.

        Args:
            name: Module name; when omitted, an ``arith_``-prefixed name is generated from
                the hash of the network payload.
            attributes: Additional module attributes, which may override the built-in keys.

        Returns:
            Operation: The public registers are all input and output ports of the network.
        """
        payload = self.payload()
        name = name or "arith_" + hashlib.sha256(payload.encode()).hexdigest()[:20]
        attrs = {
            "arithmetic_network": payload,
            "correctness": "pending",
            "workspace_contract": "zero_in_zero_out",
            **(attributes or {}),
        }
        b = Builder(
            name,
            {k: Bits(len(v)) for k, v in {**self.inputs, **self.outputs}.items()},
            attributes=attrs,
        )
        refs = {v: b[name][i] for name, bits in self.inputs.items() for i, v in enumerate(bits)}
        computed = [i for i, node in enumerate(self.nodes) if node[0] not in {"const", "input"}]
        banks = [
            b.local("ssa_" + str(i // 64), Bits(min(64, len(computed) - i)))
            for i in range(0, len(computed), 64)
        ]
        refs.update({v: banks[i // 64][i % 64] for i, v in enumerate(computed)})

        def copy(source: int, target: Ref) -> None:
            """XOR-copy the bit value of a source node into a target view; constant 1 is
            written through an X gate."""
            if source == 1:
                b.x(target)
            elif source != 0:
                b.xor(refs[source], target)

        for i in computed:
            op, *args = self.nodes[i]
            if op == "not":
                b.x(refs[i])
                copy(cast("int", args[0]), refs[i])
            elif op == "xor":
                copy(cast("int", args[0]), refs[i])
                copy(cast("int", args[1]), refs[i])
            else:
                with b.control(fuse(refs[cast("int", args[0])], refs[cast("int", args[1])])):
                    b.x(refs[i])
        forward = tuple(b._frames[0])
        for name, bits in self.outputs.items():
            for i, value in enumerate(bits):
                copy(value, b[name][i])
        from oracq.infrastructure.ir import Adjoint

        b.emit(Adjoint(forward))
        return b.finish()


DEFAULT_FIXED_FORMAT = FixedFormat()
"""Default fixed-point format: ``width=8``, ``fraction=3``, signed."""


@lru_cache(maxsize=256)
def fixed_arithmetic(kind: str, fmt: FixedFormat = DEFAULT_FIXED_FORMAT) -> Operation:
    """Outputs XOR; status[0] marks a domain failure and status[1] marks exceeding the
    word width, not a precision bound.

    Args:
        kind: Arithmetic kind, one of ``add``, ``sub``, ``neg``, ``abs``, ``mul``,
            ``div``, ``reciprocal``, ``sqrt``, ``lt``, ``eq``, ``select``, ``and``,
            ``or`` or ``xor``.
        fmt: Fixed-point format of the operands, bit width, fraction bits and signedness;
            defaults to the signed 8-bit format with 3 fraction bits.

    Returns:
        Operation: The reversible fixed-point arithmetic operation synthesized from the
        Boolean network, with output registers out and status(2).
    """
    n, f, signed = fmt.width, fmt.fraction, fmt.signed
    net = BooleanNetwork()
    a = net.input("a", n)
    unary = kind in {"neg", "abs", "sqrt", "reciprocal"}
    b = None if unary else net.input("b", n)
    invalid, overflow = 0, 0
    sign_a = a[-1] if signed else 0
    if kind == "add":
        out, carry = net.add(a, cast("list[int]", b))
        overflow = (
            net.and_(
                net.inv(net.xor(a[-1], cast("list[int]", b)[-1])), net.xor(a[-1], out[-1])
            )
            if signed
            else carry
        )
    elif kind == "sub":
        out = net.sub(a, cast("list[int]", b))
        overflow = (
            net.and_(net.xor(a[-1], cast("list[int]", b)[-1]), net.xor(a[-1], out[-1]))
            if signed
            else net.lt(a, cast("list[int]", b))
        )
    elif kind in {"neg", "abs"}:
        out = net.neg(a) if kind == "neg" else net.abs(a, signed)
        overflow = net.and_(sign_a, net.inv(net.any(a[:-1])))
    elif kind in {"mul", "div", "reciprocal"}:
        ma = net.abs(a, signed)
        if kind == "reciprocal":
            numerator, denominator, sign = net.const(1 << (2 * f), n + 2 * f + 1), ma, sign_a
        else:
            mb = net.abs(cast("list[int]", b), signed)
            sign = net.xor(sign_a, cast("list[int]", b)[-1]) if signed else 0
            numerator, denominator = [0] * f + ma, mb
        if kind == "mul":
            mag = net.mul(ma, mb, 2 * n)[f:]
        else:
            mag = net.div(numerator, denominator)
            invalid = net.inv(net.any(denominator))
            mag = net.mux(invalid, [0] * len(mag), mag)
        overflow = net.any(mag[n - int(signed) :])
        raw = net.resize(mag, n)
        out = net.mux(sign, net.neg(raw), raw)
    elif kind == "sqrt":
        invalid = sign_a
        mag = net.sqrt([0] * f + a)
        overflow = net.any(mag[n - int(signed) :])
        out = net.mux(invalid, [0] * n, net.resize(mag, n))
    elif kind in {"lt", "eq"}:
        if kind == "eq":
            bit = net.inv(
                net.any([net.xor(x, y) for x, y in zip(a, cast("list[int]", b), strict=True)])
            )
        else:
            bit = net.lt(a, cast("list[int]", b))
            if signed:
                bit = net.mux(net.xor(a[-1], cast("list[int]", b)[-1]), [a[-1]], [bit])[0]
        out = [bit]
    elif kind == "select":
        select = net.input("select", 1)[0]
        out = net.mux(select, a, cast("list[int]", b))
    elif kind in {"and", "or", "xor"}:
        fn = {"and": net.and_, "or": net.or_, "xor": net.xor}[kind]
        out = [fn(x, y) for x, y in zip(a, cast("list[int]", b), strict=True)]
    else:
        raise ValidationError("Unknown arithmetic generator: " + kind)
    net.outputs = {"out": out, "status": [invalid, overflow]}
    return net.operation(
        attributes={
            "arithmetic_kind": kind,
            "fixed_width": n,
            "fixed_fraction": f,
            "fixed_signed": signed,
            "rounding": "toward_zero; modular_wrap",
        }
    )


class BooleanCppFactory:
    """Generates a real PySparQ C++ operator for the same Boolean graph; supports
    cross-register views."""

    def __init__(self, network: BooleanNetwork, cache_dir: str) -> None:
        """Record the network, the C++ compilation cache directory, and the cache table
        from layout to operator instance."""
        self.network: BooleanNetwork = network
        self.cache_dir: str = str(cache_dir)
        self.classes: dict[
            tuple[tuple[tuple[int, int], ...], ...],
            tuple[Callable[..., object], list[tuple[str, str]]],
        ] = {}

    def __call__(self, context: NativeContext) -> object:
        """Lazily compile and instantiate the C++ Boolean operator for the register
        layout at the call site."""
        layout = tuple(
            tuple((p.start, p.width) for p in ref.parts) for ref in context.site.arguments
        )
        if layout not in self.classes:
            from pysparq.dynamic_operator import compile_operator

            net = self.network
            name = (
                "Bool_" + hashlib.sha256((net.payload() + repr(layout)).encode()).hexdigest()[:20]
            )
            fields: list[str] = []
            params: list[tuple[str, str]] = []
            assigns: list[str] = []
            read: dict[tuple[str, int], str] = {}
            write: list[str] = []
            for reg, spans in zip(context.site.module.registers, layout, strict=True):
                offset = 0
                for j, (start, width) in enumerate(spans):
                    key = "r_" + reg.name + "_" + str(j)
                    fields.append("size_t " + key + ";")
                    params.append(("size_t", key))
                    assigns.append(key + "(" + key + "_)")
                    for i in range(width):
                        read[(reg.name, offset + i)] = (
                            f"((s.registers.at({key}).value >> {start + i}) & 1ULL)"
                        )
                        if reg.name in net.outputs:
                            v = net.outputs[reg.name][offset + i]
                            write.append(
                                f"s.registers.at({key}).value ^= (uint64_t(v[{v}]) << {start + i});"
                            )
                    offset += width
            code = ["v[0]=false; v[1]=true;"]
            for i, (op, *args) in enumerate(net.nodes[2:], 2):
                expr = (
                    read[cast("tuple[str, int]", tuple(args))]
                    if op == "input"
                    else f"!v[{args[0]}]"
                    if op == "not"
                    else f"v[{args[0]}] {'^' if op == 'xor' else '&'} v[{args[1]}]"
                )
                code.append(f"v[{i}] = {expr};")
            source = (
                f"class {name} : public SelfAdjointOperator {{ "
                + " ".join(fields)
                + " public: "
                + name
                + "("
                + ", ".join(t + " " + p + "_" for t, p in params)
                + ") : "
                + ", ".join(assigns)
                + " {} "
                + "void operator()(std::vector<System>& state) const override { "
                + f"std::vector<bool> v({len(net.nodes)}); for (auto &s:state) {{ "
                + " ".join(code + write)
                + " } } };"
            )
            self.classes[layout] = (
                compile_operator(
                    name=name,
                    cpp_code=source,
                    base_class="SelfAdjointOperator",
                    constructor_args=params,
                    cache_dir=self.cache_dir,
                ),
                params,
            )
        cls, params = self.classes[layout]
        ids = [
            cast("ModuleType", context.ps).System.get_id(context.names[span.register])
            for ref in context.site.arguments
            for span in ref.parts
        ]
        return cls(**{p: rid for (_, p), rid in zip(params, ids, strict=True)})


def arithmetic_native_registry(
    program: Program, *, cache_dir: str = "out/native-cache"
) -> NativeRegistry:
    """Build a native operator registry for the arithmetic modules in a program.

    It scans modules carrying the ``arithmetic_network`` attribute, rebuilds
    the Boolean network from the attribute payload and registers a
    ``BooleanCppFactory``. C++ compilation does not happen here: the PySparQ
    operator is compiled on demand, and cached, only when the factory is
    invoked with a concrete call layout.

    Args:
        program: The program to scan.
        cache_dir: Cache directory for the C++ compilation artifacts.

    Returns:
        NativeRegistry: Registry mapping arithmetic modules to native implementations; the
        label is taken from the ``arithmetic_kind`` attribute, defaulting to
        the module name.

    Raises:
        ValidationError: The module registers and the network ports mismatch in name or
            bit width.
    """
    from oracq.infrastructure.native import NativeRegistry

    registry = NativeRegistry()
    for module in program.modules:
        payload = dict(module.attributes).get("arithmetic_network")
        if payload:
            network = BooleanNetwork.from_payload(cast("str", payload))
            ports = {r.name: r.type.width for r in module.registers}
            if ports != {k: len(v) for k, v in {**network.inputs, **network.outputs}.items()}:
                raise ValidationError("The native Boolean network does not match the module"
                                      " ports")
            registry.register(
                module,
                BooleanCppFactory(network, cache_dir),
                label=cast("str", dict(module.attributes).get("arithmetic_kind", module.name)),
            )
    return registry
