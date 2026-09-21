"""可直接调用的经典纯函数；这些函数不依赖量子 Builder。"""

from __future__ import annotations

import cmath
import math


def pressure(rho: float, momentum: float, energy: float, gamma: float = 1.4) -> float:
    """由密度、动量与总能计算理想气体压强。"""
    velocity = momentum / rho
    return (gamma - 1) * (energy - 0.5 * momentum * velocity)


def roe_speed(rho_l: float, momentum_l: float, rho_r: float, momentum_r: float) -> float:
    """计算左右状态的 Roe 平均对流速度。"""
    left = math.sqrt(rho_l)
    right = math.sqrt(rho_r)
    return (left * momentum_l / rho_l + right * momentum_r / rho_r) / (left + right)


def phase_response(z: complex) -> complex:
    """返回相位响应 exp(i z)/(1+z^2)。"""
    return cmath.exp(1j * z) / (1 + z * z)


def guarded_reciprocal(x: float) -> float:
    """带保护的倒数；输入为零时返回 0。"""
    if x == 0:
        return 0.0
    return 1.0 / x


def polynomial(x: float, order: int = 3) -> float:
    """计算 x 的 0 到 order-1 次幂之和。"""
    result = 0.0
    for k in range(order):
        result = result + x**k
    return result
