"一维 Euler frozen-Roe 矩阵元的可逆算术生成；物理/定点正确性待核验。"

from __future__ import annotations

from dataclasses import dataclass

from pyqecclang.algorithms.common.arithmetic import BooleanNetwork, FixedFormat, fixed_arithmetic
from pyqecclang.infrastructure.ir import Adjoint, Bits


@dataclass(frozen=True)
class Word:
    """定点算术表达式中的字，运算符重载后委托给所属生成器执行。

    Attributes:
        algebra: 生成该字的 ``ArithmeticBuilder``，负责实际生成算术节点。
        ref: 承载该值的寄存器引用。
    """

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
        """生成对该字求平方根的算术节点，返回承载结果的字。"""
        return self.algebra.calc("sqrt", self)

    def abs(self):
        """生成对该字求绝对值的算术节点，返回承载结果的字。"""
        return self.algebra.calc("abs", self)


class ArithmeticBuilder:
    """宿主侧表达式生成器，结束时只留下模块调用/门和私有寄存器。"""

    def __init__(self, builder, fmt):
        self.b, self.fmt = builder, fmt
        self.count, self.constants, self.flags = 0, {}, []
        self.start = len(builder._frames[0])

    def local(self, width=None):
        """分配一个新的私有寄存器，用于存放中间值或常量。

        Args:
            width: 位宽；省略时使用当前定点格式的字宽。

        Returns:
            新分配的寄存器引用，名字为 ``v_`` 加自增编号。
        """
        self.count += 1
        return self.b.local(
            "v_" + str(self.count), Bits(self.fmt.width if width is None else width)
        )

    def word(self, value):
        """把各类输入统一包装成 ``Word``。

        ``Word`` 原样返回；带 ``parts`` 属性的寄存器引用直接包装；数值按当前
        定点格式编码，相同编码的常量只物化一次并复用同一寄存器。

        Args:
            value: ``Word``、寄存器引用或可按定点格式编码的数值。

        Returns:
            Word: 与输入等价的算术字。
        """
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
        """生成一个定点算术运算节点并返回其输出字。

        同一寄存器被重复用作操作数时，先用 XOR 复制到新寄存器、调用结束后
        再还原，以支持平方等自身参与运算的情形；节点的状态寄存器记入
        ``flags``，由 ``finish`` 统一汇总。

        Args:
            kind: 运算类型名，如 ``add``、``mul``、``sqrt``。
            *args: 操作数，可为 ``Word``、寄存器或数值。
            select: ``select`` 运算的控制位，其他运算不使用。

        Returns:
            Word: 承载运算结果的字。
        """
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
        """按单个控制位在两个字之间选择，位值 1 取 ``yes``、0 取 ``no``。"""
        return self.calc("select", yes, no, select=bit)

    def choose3(self, index, values):
        """用两位索引在三个候选字中选择。

        Args:
            index: 两位控制位序列，低位在前。
            values: 三个候选字，依次对应索引 0..2。

        Returns:
            Word: 被选中的字；索引为 3 时结果为零。
        """
        low = self.choose(index[0], values[1], values[0])
        high = self.choose(index[0], 0, values[2])
        return self.choose(index[1], high, low)

    def finish(self, outputs, status):
        """收尾并返回最终操作。

        每个输出字 XOR 进对应目标寄存器；各算术节点的状态位按位 OR 写入
        ``status``；构造期间发出的正向指令以 ``Adjoint`` 反向重放，消除全部
        中间值。

        Args:
            outputs: 目标寄存器与结果字的配对序列。
            status: 接收聚合状态的两位寄存器。

        Returns:
            Operation: 生成的完整操作。
        """
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

    from pyqecclang.applications.roe_formulas import frozen_roe_face
    from pyqecclang.infrastructure.mathfunc import Index, compile_function

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
