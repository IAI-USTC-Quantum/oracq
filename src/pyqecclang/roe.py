"""一维 Euler frozen-Roe 矩阵元的可逆算术生成；物理/定点正确性待核验。"""

from __future__ import annotations

from dataclasses import dataclass

from .arithmetic import BooleanNetwork, FixedFormat, fixed_arithmetic
from .ir import Adjoint, Bits


@dataclass(frozen=True)
class Word:
    algebra: object
    ref: object

    def __add__(self, other):
        return self.algebra.calc("add", self, other)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self.algebra.calc("sub", self, other)

    def __rsub__(self, other):
        return self.algebra.word(other) - self

    def __mul__(self, other):
        return self.algebra.calc("mul", self, other)

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        return self.algebra.calc("div", self, other)

    def __neg__(self):
        return self.algebra.calc("neg", self)

    def sqrt(self):
        return self.algebra.calc("sqrt", self)

    def abs(self):
        return self.algebra.calc("abs", self)


class ArithmeticBuilder:
    """宿主侧表达式生成器，结束时只留下模块调用/门和私有寄存器。"""

    def __init__(self, builder, fmt):
        self.b, self.fmt = builder, fmt
        self.count, self.constants, self.flags = 0, {}, []
        self.start = len(builder._frames[0])

    def local(self, width=None):
        self.count += 1
        return self.b.local(
            "v_" + str(self.count), Bits(self.fmt.width if width is None else width)
        )

    def word(self, value):
        if isinstance(value, Word):
            return value
        if hasattr(value, "parts"):
            return Word(self, value)
        raw = self.fmt.encode(value)
        if raw not in self.constants:
            ref = self.local()
            for i in range(ref.width):
                if (raw >> i) & 1:
                    self.b.x(ref[i])
            self.constants[raw] = Word(self, ref)
        return self.constants[raw]

    def calc(self, kind, *args, select=None):
        op = fixed_arithmetic(kind, self.fmt)
        refs, copied = [], []
        for arg in args:
            ref = self.word(arg).ref
            if ref in refs:
                clone = self.local()
                self.b.xor(ref, clone)
                copied.append((ref, clone))
                ref = clone
            refs.append(ref)
        out_width = next(r.type.width for r in op.module.registers if r.name == "out")
        out, status = self.local(out_width), self.local(2)
        kwargs = dict(zip(("a", "b"), refs, strict=False))
        kwargs.update(out=out, status=status)
        if select is not None:
            kwargs["select"] = select
        self.b.call(op, **kwargs)
        for source, clone in reversed(copied):
            self.b.xor(source, clone)
        self.flags.append(status)
        return Word(self, out)

    def choose(self, bit, yes, no):
        return self.calc("select", yes, no, select=bit)

    def choose3(self, index, values):
        low = self.choose(index[0], values[1], values[0])
        high = self.choose(index[0], 0, values[2])
        return self.choose(index[1], high, low)

    def finish(self, outputs, status):
        outputs = [(target, self.word(word)) for target, word in outputs]
        forward = tuple(self.b._frames[0][self.start :])
        for target, word in outputs:
            self.b.xor(word.ref, target)
        # 输出状态为各算术节点状态的 OR，不以 XOR 抵消重复错误。
        if self.flags:
            net = BooleanNetwork()
            all_flags = [net.input("f_" + str(i), 2) for i in range(len(self.flags))]
            net.outputs = {"status": [net.any([f[j] for f in all_flags]) for j in range(2)]}
            op = net.operation()
            self.b.call(
                op, status=status, **{"f_" + str(i): flag for i, flag in enumerate(self.flags)}
            )
        self.b.emit(Adjoint(forward))
        return self.b.finish()


def roe_face(*, fmt=None, gamma=1.4, entropy_delta=0.125):
    """普通 Python Roe 函数自动编译，保留现有六个场量及行列输入 ABI。"""
    from dataclasses import replace

    from .mathfunc import Index, compile_function
    from .mathfunc.roe_formulas import frozen_roe_face

    fmt = fmt or FixedFormat(10, 5)
    compiled = compile_function(
        frozen_roe_face,
        fmt=fmt,
        inputs={
            **{key: "real" for key in ("rho_l", "m_l", "e_l", "rho_r", "m_r", "e_r")},
            "row": Index(2),
            "col": Index(2),
        },
        constants={"gamma": gamma, "entropy_delta": entropy_delta},
        output_names=("left", "right"),
    )
    attrs = dict(compiled.operation.module.attributes)
    attrs.update(
        algorithm="frozen_roe_1d_euler",
        matrix_entry_source="compiled_pure_function_on_raw_conserved_variables",
    )
    return replace(
        compiled.operation,
        module=replace(compiled.operation.module, attributes=tuple(sorted(attrs.items()))),
    )
