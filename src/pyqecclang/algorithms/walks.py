"""离散时间 coined quantum walk。"""

from pyqecclang.algorithms.contracts import positive_integer
from pyqecclang.algorithms.operators import _name
from pyqecclang.infrastructure.builder import Builder
from pyqecclang.infrastructure.ir import Bits


def cycle_walk(width, *, steps=1):
    """生成周期格点上的 Hadamard coined walk。

    Args:
        width: position 位宽，周期长度为 2**width。
        steps: 非负步数。

    Returns:
        Operation: 公开 position 和一位 coin。每步先更新 coin，再按 0/1 条件分别移动 +1/-1。

    输入态由调用方准备；零输入对应从位置零、coin 零开始。"""
    positive_integer(width, "cycle_walk.width", maximum=64)
    positive_integer(steps, "cycle_walk.steps", minimum=0)
    b = Builder(
        _name("cycle_walk", width, steps),
        {"position": Bits(width), "coin": Bits(1)},
        attributes={"algorithm": "coined_cycle_walk", "steps": steps},
    )
    with b.repeat(steps):
        b.h(b["coin"])
        with b.control(b["coin"], 0):
            b.add_const(b["position"].reinterpret("uint"), 1)
        with b.control(b["coin"], 1):
            b.add_const(b["position"].reinterpret("uint"), (1 << width) - 1)
    return b.finish()
