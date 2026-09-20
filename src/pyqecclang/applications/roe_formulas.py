"普通、无副作用的 Roe 数学函数；可经典调用，也可由 compile_function 编译。"

import math


def pick3(index, first, second, third):
    """按整数索引选取三个值之一，索引不在 0..2 范围时返回 ``0.0``。"""
    if index == 0:
        return first
    if index == 1:
        return second
    if index == 2:
        return third
    return 0.0


def conserved_to_primitive(rho, momentum, energy, gamma):
    """把一维 Euler 守恒变量换算为原始变量。

    Args:
        rho: 密度。
        momentum: 动量密度。
        energy: 总能量密度。
        gamma: 比热比。

    Returns:
        tuple: ``(velocity, pressure, enthalpy)``，焓为总能量加压强除以密度。
    """
    velocity = momentum / rho
    pressure = (gamma - 1) * (energy - 0.5 * momentum * velocity)
    enthalpy = (energy + pressure) / rho
    return velocity, pressure, enthalpy


def euler_entry(velocity, enthalpy, row, col, gamma):
    """返回一维 Euler 通量 Jacobi 矩阵的指定元素。

    元素只用流速与焓表示；行、列索引与守恒变量 ``(rho, momentum, energy)``
    的次序对应，从零开始编号。

    Args:
        velocity: 流速。
        enthalpy: 比焓，定义同 ``conserved_to_primitive`` 的返回值。
        row: 行索引，取 0..2。
        col: 列索引，取 0..2。
        gamma: 比热比。

    Returns:
        float: Jacobi 矩阵第 ``row`` 行、第 ``col`` 列的元素。
    """
    square = velocity * velocity
    first = pick3(col, 0.0, 1.0, 0.0)
    second = pick3(col, 0.5 * (gamma - 3) * square, (3 - gamma) * velocity, gamma - 1)
    third = pick3(
        col,
        velocity * (0.5 * (gamma - 1) * square - enthalpy),
        enthalpy - (gamma - 1) * square,
        gamma * velocity,
    )
    return pick3(row, first, second, third)


def entropy_absolute(eigenvalue, delta):
    """带 Harten 熵修正的特征值绝对值。

    ``delta`` 为正且特征值绝对值小于它时，返回平滑值
    ``(eigenvalue**2+delta**2)/(2*delta)``；否则直接返回绝对值。

    Args:
        eigenvalue: 特征值。
        delta: 熵修正阈值，非正时关闭修正。

    Returns:
        float: 修正后的特征值绝对值。
    """
    value = abs(eigenvalue)
    if delta <= 0:
        return value
    if value >= delta:
        return value
    return (eigenvalue * eigenvalue + delta * delta) / (2 * delta)


def frozen_roe_face(rho_l, m_l, e_l, rho_r, m_r, e_r, row, col, gamma=1.4, entropy_delta=0.125):
    """计算 frozen-Roe 左右系数矩阵在指定行列处的元素。

    由左右守恒状态做 Roe 平均得到平均流速与焓；特征值 ``u-c``、``u``、
    ``u+c`` 经熵修正取绝对值后，按 ``R*diag(a)*R**-1`` 组装 ``|A|`` 的元素。
    再与 ``euler_entry`` 给出的两侧通量 Jacobi 矩阵元素组合，返回
    ``0.5*(A_l+|A|)`` 与 ``0.5*(A_r-|A|)`` 在 ``row``、``col`` 处的值，作为
    Roe 通量对左右状态的分摊系数。

    Args:
        rho_l: 左侧密度。
        m_l: 左侧动量密度。
        e_l: 左侧总能量密度。
        rho_r: 右侧密度。
        m_r: 右侧动量密度。
        e_r: 右侧总能量密度。
        row: 矩阵元素行索引，取 0..2。
        col: 矩阵元素列索引，取 0..2。
        gamma: 比热比。
        entropy_delta: Harten 熵修正阈值。

    Returns:
        tuple: ``(left, right)``，左右系数矩阵在指定行列处的元素。
    """
    ul, pl, hl = conserved_to_primitive(rho_l, m_l, e_l, gamma)
    ur, pr, hr = conserved_to_primitive(rho_r, m_r, e_r, gamma)
    wl = math.sqrt(rho_l)
    wr = math.sqrt(rho_r)
    total = wl + wr
    u = (wl * ul + wr * ur) / total
    h = (wl * hl + wr * hr) / total
    u2 = u * u
    c2 = (gamma - 1) * (h - 0.5 * u2)
    c = math.sqrt(c2)
    beta = (gamma - 1) / c2
    inverse_c = 0.5 / c
    bu = beta * u
    bu2 = beta * u2
    r0 = pick3(row, 1.0, u - c, h - u * c)
    r1 = pick3(row, 1.0, u, 0.5 * u2)
    r2 = pick3(row, 1.0, u + c, h + u * c)
    inv0 = pick3(col, 0.25 * bu2 + u * inverse_c, -0.5 * bu - inverse_c, 0.5 * beta)
    inv1 = pick3(col, 1 - 0.5 * bu2, bu, -beta)
    inv2 = pick3(col, 0.25 * bu2 - u * inverse_c, -0.5 * bu + inverse_c, 0.5 * beta)
    a0 = entropy_absolute(u - c, entropy_delta)
    a1 = entropy_absolute(u, entropy_delta)
    a2 = entropy_absolute(u + c, entropy_delta)
    absolute = r0 * a0 * inv0 + r1 * a1 * inv1 + r2 * a2 * inv2
    left = 0.5 * (euler_entry(ul, hl, row, col, gamma) + absolute)
    right = 0.5 * (euler_entry(ur, hr, row, col, gamma) - absolute)
    return left, right
