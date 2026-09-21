"可逆定点算术：Boolean SSA → compute/XOR/uncompute，及同图原生实现。"

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache

from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits, ValidationError, fuse


@dataclass(frozen=True)
class FixedFormat:
    """描述可逆定点算术使用的定点数格式。

    ``width`` 位字长的低 ``fraction`` 位为小数位；``signed`` 为 True 时最高位
    为补码符号位。数值与位模式之间通过 ``encode`` 与 ``decode`` 换算。

    Attributes:
        width: 总位宽，含符号位（若存在），有效范围为 2..64。
        fraction: 小数位数，不得占用符号位。
        signed: 是否以最高位为补码符号位。

    Raises:
        ValidationError: ``width`` 越出 2..64，或 ``fraction`` 占用到符号位。
    """
    width: int = 8
    fraction: int = 3
    signed: bool = True

    def __post_init__(self):
        if not 2 <= self.width <= 64 or not 0 <= self.fraction < self.width - int(self.signed):
            raise ValidationError("定点格式要求 2≤width≤64 且小数位不占用符号位")

    def encode(self, value):
        """把定点数值编码为 ``width`` 位整数。

        数值先乘以 ``2**fraction`` 并向零截断，再掩码到 ``width`` 位。

        Args:
            value: 待编码的定点数值。

        Returns:
            int: 按 ``2**width`` 回绕后的整数位模式。
        """
        return int(value * (1 << self.fraction)) & ((1 << self.width) - 1)

    def decode(self, value):
        """把 ``width`` 位整数位模式解码回定点数值。

        Args:
            value: 整数位模式，越出字长的高位先被掩掉。

        Returns:
            float: 有符号格式按补码解释符号位后，再除以 ``2**fraction`` 的值。
        """
        value &= (1 << self.width) - 1
        if self.signed and value >> (self.width - 1):
            value -= 1 << self.width
        return value / (1 << self.fraction)


class BooleanNetwork:
    """0/1 为常量；其余编号是输入或不可变布尔值。"""

    def __init__(self):
        self.nodes = [("const", 0), ("const", 1)]
        self.cache = {}
        self.inputs = {}
        self.outputs = {}

    def node(self, key):
        """返回键对应的节点编号，不存在时创建新节点。

        以元组键做 hash-consing：相同键始终返回相同编号。

        Args:
            key: 描述常量、输入或运算的可迭代键。

        Returns:
            int: 该键的节点编号。
        """
        key = tuple(key)
        if key not in self.cache:
            self.cache[key] = len(self.nodes)
            self.nodes.append(key)
        return self.cache[key]

    def input(self, name, width):
        """声明一个输入端口并登记其逐位节点。

        Args:
            name: 端口名。
            width: 端口位宽。

        Returns:
            list: 节点编号列表，下标为位序（低位在前）。
        """
        bits = [self.node(("input", name, i)) for i in range(width)]
        self.inputs[name] = bits
        return bits

    def inv(self, a):
        """对节点取反，带常量折叠与双重取消化简。

        Args:
            a: 节点编号。

        Returns:
            int: 取反结果的节点编号。
        """
        if a < 2:
            return 1 - a
        if self.nodes[a][0] == "not":
            return self.nodes[a][1]
        return self.node(("not", a))

    def xor(self, a, b):
        """异或两个节点，带常量单位元与 ``a xor a = 0`` 化简。

        Args:
            a: 节点编号。
            b: 节点编号。

        Returns:
            int: 异或结果的节点编号。
        """
        if a == b:
            return 0
        if a == 0 or b == 0:
            return a or b
        if a == 1 or b == 1:
            return self.inv(b if a == 1 else a)
        return self.node(("xor", *sorted((a, b))))

    def and_(self, a, b):
        """对两个节点做按位与，带常量吸收与幂等化简。

        Args:
            a: 节点编号。
            b: 节点编号。

        Returns:
            int: 与结果的节点编号。
        """
        if a == b:
            return a
        if a == 0 or b == 0:
            return 0
        if a == 1 or b == 1:
            return b if a == 1 else a
        return self.node(("and", *sorted((a, b))))

    def or_(self, a, b):
        """对两个节点做按位或，经 ``xor`` 与 ``and_`` 复合实现。

        Args:
            a: 节点编号。
            b: 节点编号。

        Returns:
            int: 或结果的节点编号。
        """
        return self.xor(self.xor(a, b), self.and_(a, b))

    def any(self, bits):
        """把一组节点归约为它们的按位或。

        Args:
            bits: 节点编号的可迭代对象。

        Returns:
            int: 归约结果的节点编号；空输入返回常量 0。
        """
        result = 0
        for bit in bits:
            result = self.or_(result, bit)
        return result

    def mux(self, select, yes, no):
        """按选择位在两组等宽信号间逐位选择。

        Args:
            select: 选择位节点编号。
            yes: 选择位为 1 时输出的位列表。
            no: 选择位为 0 时输出的位列表，须与 ``yes`` 同宽。

        Returns:
            list: 逐位选择结果的节点编号列表。
        """
        return [
            self.xor(b, self.and_(select, self.xor(a, b))) for a, b in zip(yes, no, strict=True)
        ]

    @staticmethod
    def const(value, width):
        """把非负整数展开为 ``width`` 位常量位列表。

        Args:
            value: 非负整数常量。
            width: 位宽。

        Returns:
            list: 各位取 0 或 1，低位在前。
        """
        return [(value >> i) & 1 for i in range(width)]

    @staticmethod
    def resize(bits, width):
        """把位列表截断或高位补零到指定位宽。

        Args:
            bits: 原位列表。
            width: 目标位宽。

        Returns:
            list: 长度恰为 ``width`` 的新列表。
        """
        return (list(bits) + [0] * width)[:width]

    def add(self, a, b, carry=0):
        """行波进位加法。

        Args:
            a: 加数的位列表。
            b: 与 ``a`` 同宽的加数位列表。
            carry: 进入最低位的初始进位节点编号。

        Returns:
            tuple: ``(result, carry)``，逐位和与最高位向前的进位。
        """
        result = []
        for x, y in zip(a, b, strict=True):
            p = self.xor(x, y)
            result.append(self.xor(p, carry))
            carry = self.xor(self.and_(x, y), self.and_(p, carry))
        return result, carry

    def neg(self, a):
        """补码取负：按位取反后在最低位进一。

        Args:
            a: 操作数的位列表。

        Returns:
            list: ``-a`` 的位列表。
        """
        return self.add([self.inv(x) for x in a], [0] * len(a), 1)[0]

    def sub(self, a, b):
        """减法，按 ``a + ~b + 1`` 实现。

        Args:
            a: 被减数的位列表。
            b: 与 ``a`` 同宽的减数位列表。

        Returns:
            list: ``a - b`` 的位列表。
        """
        return self.add(a, [self.inv(x) for x in b], 1)[0]

    def lt(self, a, b):
        """无符号比较 ``a < b``。

        Args:
            a: 左操作数的位列表。
            b: 与 ``a`` 同宽的右操作数位列表。

        Returns:
            int: 比较结果的节点编号，真值为 1。
        """
        return self.inv(self.add(a, [self.inv(x) for x in b], 1)[1])

    def abs(self, a, signed):
        """返回操作数的绝对值。

        Args:
            a: 操作数的位列表。
            signed: 为 True 时按符号位在 ``a`` 与 ``-a`` 间选择；为 False 时原样返回。

        Returns:
            list: 绝对值的位列表。
        """
        return self.mux(a[-1], self.neg(a), a) if signed else a

    def mul(self, a, b, width):
        """移位累加乘法，乘积截断到 ``width`` 位。

        Args:
            a: 被乘数的位列表。
            b: 乘数的位列表。
            width: 结果位宽。

        Returns:
            list: 乘积低 ``width`` 位的位列表。
        """
        result = [0] * width
        for i, bit in enumerate(b[:width]):
            row = [0] * i + [self.and_(bit, x) for x in a[: width - i]]
            result = self.add(result, self.resize(row, width))[0]
        return result

    def div(self, a, b):
        """恢复余数法除法，只返回商。

        被除数从最高位起逐位移入部分余数；余数不小于除数时减去除数并置商位。

        Args:
            a: 被除数的位列表。
            b: 除数的位列表。

        Returns:
            list: 商的位列表，位宽与 ``a`` 相同；余数被丢弃。
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

    def sqrt(self, a):
        """逐位恢复式开方。

        输入位数为奇数时先高位补零到偶数，被开方数按两位一组从高位移入。

        Args:
            a: 被开方数的位列表。

        Returns:
            list: 平方根的位列表，位数为补零后输入位数的一半。
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

    def evaluate(self, **inputs):
        """在经典侧对整张网络求值。

        Args:
            **inputs: 端口名到非负整数的映射，按位读入各输入端口。

        Returns:
            dict: 每个输出端口名到其整数求值结果。
        """
        values = [0, 1]
        for op, *args in self.nodes[2:]:
            if op == "input":
                name, bit = args
                v = (inputs[name] >> bit) & 1
            elif op == "not":
                v = 1 ^ values[args[0]]
            elif op == "xor":
                v = values[args[0]] ^ values[args[1]]
            else:
                v = values[args[0]] & values[args[1]]
            values.append(v)
        return {
            name: sum(values[v] << i for i, v in enumerate(bits))
            for name, bits in self.outputs.items()
        }

    def payload(self):
        """把网络序列化为紧凑 JSON 字符串。

        Returns:
            str: 含节点、输入与输出端口的载荷，可作为缓存键或模块属性。
        """
        return json.dumps(
            dict(nodes=self.nodes, inputs=self.inputs, outputs=self.outputs), separators=(",", ":")
        )

    @classmethod
    def from_payload(cls, value):
        """从 JSON 载荷重建网络，并把载荷当不可信输入做完整校验。

        校验覆盖常量前缀、端口命名与互斥、位宽 1..64、位引用范围、
        顺序有向无环性质以及输入端口映射的一致性。

        Args:
            value: ``payload`` 生成的 JSON 字符串。

        Returns:
            BooleanNetwork: 校验通过后重建的网络。

        Raises:
            ValidationError: 载荷任一结构或一致性检查失败。
        """
        data = json.loads(value)
        net = cls()
        net.nodes = [tuple(x) for x in data["nodes"]]
        net.inputs, net.outputs = data["inputs"], data["outputs"]
        from pyqecclang.infrastructure.validation import name

        if net.nodes[:2] != [("const", 0), ("const", 1)] or set(net.inputs) & set(net.outputs):
            raise ValidationError("非法 Boolean 网络常量或端口")
        for port, bits in {**net.inputs, **net.outputs}.items():
            name(port)
            if not isinstance(bits, list) or not 1 <= len(bits) <= 64:
                raise ValidationError("Boolean 网络端口位宽无效")
            if any(type(v) is not int or not 0 <= v < len(net.nodes) for v in bits):
                raise ValidationError("Boolean 网络位引用越界")
        for index, node in enumerate(net.nodes[2:], 2):
            op, *args = node
            if op == "input":
                if len(args) != 2 or args[0] not in net.inputs or type(args[1]) is not int:
                    raise ValidationError("Boolean 网络输入节点无效")
                if (
                    not 0 <= args[1] < len(net.inputs[args[0]])
                    or net.inputs[args[0]][args[1]] != index
                ):
                    raise ValidationError("Boolean 网络输入映射不一致")
            elif op not in {"not", "and", "xor"} or len(args) != (1 if op == "not" else 2):
                raise ValidationError("未知 Boolean 网络节点")
            elif any(type(v) is not int or not 0 <= v < index for v in args):
                raise ValidationError("Boolean 网络不是顺序有向无环图")
        for port, bits in net.inputs.items():
            if any(net.nodes[v] != ("input", port, i) for i, v in enumerate(bits)):
                raise ValidationError("Boolean 网络输入端口不完整")
        return net

    def operation(self, name=None, *, attributes=None):
        """把网络编译为 compute/copy/uncompute 形式的可逆量子操作。

        每个运算节点先正向写入按 64 位分组的 ``ssa_`` 前缀私有 bank
        （``not`` 与 ``xor`` 用 X/XOR 拷贝实现，``and`` 用受控 X 实现），
        输出位 XOR 拷贝到公开寄存器后，再发出正向计算帧的伴随以复净全部
        bank，满足 ``zero_in_zero_out`` 工作区契约。

        Args:
            name: 模块名；省略时由网络载荷哈希生成的 ``arith_`` 前缀名。
            attributes: 追加的模块属性，可覆盖内置键。

        Returns:
            Operation: 公开寄存器为网络的全部输入与输出端口。
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

        def copy(source, target):
            if source == 1:
                b.x(target)
            elif source != 0:
                b.xor(refs[source], target)

        for i in computed:
            op, *args = self.nodes[i]
            if op == "not":
                b.x(refs[i])
                copy(args[0], refs[i])
            elif op == "xor":
                copy(args[0], refs[i])
                copy(args[1], refs[i])
            else:
                with b.control(fuse(refs[args[0]], refs[args[1]])):
                    b.x(refs[i])
        forward = tuple(b._frames[0])
        for name, bits in self.outputs.items():
            for i, value in enumerate(bits):
                copy(value, b[name][i])
        from pyqecclang.infrastructure.ir import Adjoint

        b.emit(Adjoint(forward))
        return b.finish()


DEFAULT_FIXED_FORMAT = FixedFormat()
"""默认定点格式：``width=8``、``fraction=3``、有符号。"""


@lru_cache(maxsize=256)
def fixed_arithmetic(kind, fmt=DEFAULT_FIXED_FORMAT):
    """输出 XOR；status[0]=定义域失效，status[1]=越出字长（非精度界）。"""
    n, f, signed = fmt.width, fmt.fraction, fmt.signed
    net = BooleanNetwork()
    a = net.input("a", n)
    unary = kind in {"neg", "abs", "sqrt", "reciprocal"}
    b = None if unary else net.input("b", n)
    invalid, overflow = 0, 0
    sign_a = a[-1] if signed else 0
    if kind == "add":
        out, carry = net.add(a, b)
        overflow = (
            net.and_(net.inv(net.xor(a[-1], b[-1])), net.xor(a[-1], out[-1])) if signed else carry
        )
    elif kind == "sub":
        out = net.sub(a, b)
        overflow = (
            net.and_(net.xor(a[-1], b[-1]), net.xor(a[-1], out[-1])) if signed else net.lt(a, b)
        )
    elif kind in {"neg", "abs"}:
        out = net.neg(a) if kind == "neg" else net.abs(a, signed)
        overflow = net.and_(sign_a, net.inv(net.any(a[:-1])))
    elif kind in {"mul", "div", "reciprocal"}:
        ma = net.abs(a, signed)
        if kind == "reciprocal":
            numerator, denominator, sign = net.const(1 << (2 * f), n + 2 * f + 1), ma, sign_a
        else:
            mb = net.abs(b, signed)
            sign = net.xor(sign_a, b[-1]) if signed else 0
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
            bit = net.inv(net.any([net.xor(x, y) for x, y in zip(a, b, strict=True)]))
        else:
            bit = net.lt(a, b)
            if signed:
                bit = net.mux(net.xor(a[-1], b[-1]), [a[-1]], [bit])[0]
        out = [bit]
    elif kind == "select":
        select = net.input("select", 1)[0]
        out = net.mux(select, a, b)
    elif kind in {"and", "or", "xor"}:
        fn = {"and": net.and_, "or": net.or_, "xor": net.xor}[kind]
        out = [fn(x, y) for x, y in zip(a, b, strict=True)]
    else:
        raise ValidationError("未知算术生成器：" + kind)
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
    """为同一个 Boolean 图生成真实 PySparQ C++ 算子；支持跨寄存器视图。"""

    def __init__(self, network, cache_dir):
        self.network, self.cache_dir, self.classes = network, str(cache_dir), {}

    def __call__(self, context):
        layout = tuple(
            tuple((p.start, p.width) for p in ref.parts) for ref in context.site.arguments
        )
        if layout not in self.classes:
            from pysparq.dynamic_operator import compile_operator

            net = self.network
            name = (
                "Bool_" + hashlib.sha256((net.payload() + repr(layout)).encode()).hexdigest()[:20]
            )
            fields, params, assigns, read, write = [], [], [], {}, []
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
                    read[tuple(args)]
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
            context.ps.System.get_id(context.names[span.register])
            for ref in context.site.arguments
            for span in ref.parts
        ]
        return cls(**{p: rid for (_, p), rid in zip(params, ids, strict=True)})


def arithmetic_native_registry(program, *, cache_dir="out/native-cache"):
    """为程序中的算术模块构建原生算子注册表。

    扫描带 ``arithmetic_network`` 属性的模块，按属性载荷重建 Boolean 网络
    并注册 ``BooleanCppFactory``。C++ 编译不在此处发生：PySparQ 算子在工厂
    随具体调用布局被调用时才按需编译并缓存。

    Args:
        program: 待扫描的程序。
        cache_dir: C++ 编译产物的缓存目录。

    Returns:
        NativeRegistry: 算术模块到原生实现的注册表；标签取 ``arithmetic_kind``
        属性，缺省为模块名。

    Raises:
        ValidationError: 模块寄存器与网络端口的名称或位宽不一致。
    """
    from pyqecclang.infrastructure.native import NativeRegistry

    registry = NativeRegistry()
    for module in program.modules:
        payload = dict(module.attributes).get("arithmetic_network")
        if payload:
            network = BooleanNetwork.from_payload(payload)
            ports = {r.name: r.type.width for r in module.registers}
            if ports != {k: len(v) for k, v in {**network.inputs, **network.outputs}.items()}:
                raise ValidationError("原生 Boolean 网络与模块端口不匹配")
            registry.register(
                module,
                BooleanCppFactory(network, cache_dir),
                label=dict(module.attributes).get("arithmetic_kind", module.name),
            )
    return registry
